/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/
// SPDX-License-Identifier: GPL-3.0-or-later

#include <atomic>
#include <numeric>
#include <thread>
#include "masking.h"
#include "tantan/LambdaCalculator.hh"
#include "tantan.h"
#include "blast/blast_filter.h"
#include "data/sequence_set.h"
#include "basic/config.h"
#include "util/math/math.h"

using std::unique_ptr;
using std::atomic;
using std::thread;
using std::vector;
using std::pair;
using std::atomic_size_t;

const EMap<MaskingAlgo> EnumTraits<MaskingAlgo>::to_string{ {MaskingAlgo::NONE, "None"}, {MaskingAlgo::SEG, "SEG"}, {MaskingAlgo::TANTAN, "tantan"} };
const SEMap<MaskingAlgo> EnumTraits<MaskingAlgo>::from_string{
	{ "0", MaskingAlgo::NONE },
	{ "none", MaskingAlgo::NONE },
	{ "seg", MaskingAlgo::SEG },
	{ "tantan", MaskingAlgo::TANTAN }
};
const SEMap<MaskingMode> EnumTraits<MaskingMode>::from_string{
	{"0", {MaskingMode::NONE, false}},
	{"none", MaskingMode::NONE},
	{"1", {MaskingMode::TANTAN, false}},
	{"tantan", MaskingMode::TANTAN},
	{"seg", MaskingMode::BLAST_SEG},
	{"seg-all", MaskingMode::BLAST_SEG_ALL }
};

unique_ptr<Masking> Masking::instance;

MaskingTable::MaskingTable():
	seq_count_(0),
	masked_letters_(0)
{}

MaskingTable& MaskingTable::operator=(const MaskingTable& t) {
	seq_count_ = t.seq_count_;
	masked_letters_ = t.masked_letters_;
	entry_ = t.entry_;
	seqs_ = t.seqs_;
	return *this;
}

MaskingTable::MaskingTable(const MaskingTable& t):
	seq_count_(t.seq_count_),
	masked_letters_(t.masked_letters_),
	entry_(t.entry_),
	seqs_(t.seqs_)
{}

bool MaskingTable::blank() const {
	return seq_count_ == 0;
}

size_t MaskingTable::masked_letters() const {
	return masked_letters_;
}

void MaskingTable::add(const size_t block_id, const Loc begin, const Loc end, Letter* seq) {
	{
		std::lock_guard<std::mutex> lock(mtx_);
		entry_.emplace_back(block_id, begin);
		seqs_.push_back(seq + begin, seq + end);
		++seq_count_;
		masked_letters_ += end - begin;
	}
	std::fill(seq + begin, seq + end, MASK_LETTER);
}

void MaskingTable::remove(SequenceSet& seqs, const int template_len, const bool add_bit_mask) const {
	for (size_t i = 0; i < entry_.size(); ++i) {
		Letter* ptr = seqs.ptr(entry_[i].block_id);
		std::copy(seqs_.ptr(i), seqs_.end(i), ptr + entry_[i].begin);
		if (add_bit_mask) {
			const int i0 = std::max(entry_[i].begin - template_len + 1, 0), i1 = entry_[i].begin + (int)seqs_.length(i);
			for (int j = i0; j < i1; ++j)
				ptr[j] |= SEED_MASK;
		}
	}
}

void MaskingTable::apply(SequenceSet& seqs) const {
	for (size_t i = 0; i < entry_.size(); ++i) {
		Letter* ptr = seqs.ptr(entry_[i].block_id) + entry_[i].begin;
		std::fill(ptr, ptr + seqs_.length(i), MASK_LETTER);
	}
}

static size_t mask_motifs(Letter* seq, const size_t len, const size_t block_id, MaskingTable& table) {
	if (len < MOTIF_LEN)
		return 0;
	Mask::Ranges pos;
	KmerIterator<MOTIF_LEN> it(Sequence(seq, len));
	while (it.good()) {
		if (motif_table.find(*it) != motif_table.end())
			pos.push_back(it - seq, it - seq + MOTIF_LEN);
		++it;
	}
	const ptrdiff_t n = std::accumulate(pos.cbegin(), pos.cend(), (ptrdiff_t)0, [](const ptrdiff_t s, const pair<ptrdiff_t, ptrdiff_t>& r) { return s + r.second - r.first; });

	if ((double)n / len >= 0.5)
		return 0;

	for (auto i = pos.cbegin(); i != pos.cend(); ++i)
		if (i->second - i->first <= config.max_motif_len)
			table.add(block_id, i->first, i->second, seq);

	return n;
}

Masking::Masking(const ScoreMatrix& score_matrix)
{
	const unsigned n = value_traits.alphabet_size;
	int int_matrix[20][20], *int_matrix_ptr[20];
	std::copy(int_matrix, int_matrix + 20, int_matrix_ptr);
	for (size_t i = 0; i < 20; ++i)
		for (size_t j = 0; j < 20; ++j)
			int_matrix[i][j] = score_matrix(i, j);
	cbrc::LambdaCalculator lc;
	lc.calculate(int_matrix_ptr, 20);
	
	const double lambda = lc.lambda(); // 0.324032
	for (size_t i = 0; i < size; ++i) {
		mask_table_x_[i] = value_traits.mask_char;
		mask_table_bit_[i] = (int8_t)i | bit_mask;
		for (size_t j = 0; j < size; ++j)
			if (i < n && j < n) {
				likelihoodRatioMatrixf_[i][j] = (float)Math::exp(lambda * score_matrix(i, j));
			}
	}
	std::copy(likelihoodRatioMatrixf_, likelihoodRatioMatrixf_ + size, probMatrixPointersf_);

	blast_seg_ = SegParametersNewAa();
}

Masking::~Masking() {
	SegParametersFree(blast_seg_);
}

extern double lnfact[];
namespace {
constexpr int kSegDelta[10] = {0, 200000, 275489, 324511, 360964, 390014, 414170, 434852, 452933, 468995};
constexpr int kSegLoCut[3] = {1521928, 1232933, 960000}, kSegHiCut[3] = {1221928, 962933, 720000};
inline void seg_dec_sv(int32_t* sv, int32_t cl) { for (int32_t s; (s = *sv++) != 0;) if (s == cl && *sv < cl) { sv[-1] = s - 1; break; } }
inline void seg_inc_sv(int32_t* sv, int32_t cl) { for (;;) if (*sv++ == cl) { sv[-1]++; break; } }
inline double seg_lnfact(uint32_t n) { return n < 10001u ? lnfact[n] : ((n + 0.5) * Math::log(n) - n + 0.9189385332); }
double seg_prob(const int32_t* sv, int32_t total) {
	const double totseq = ((double)total) * 2.9957322735539909;
	double ans1 = lnfact[20], ans2 = 0.0;
	if (sv[0] != 0) {
		int rem = 20, cl = 1, svim1 = sv[0], svi = sv[0];
		for (const int32_t* p = sv;; svim1 = svi) {
			if (++p - sv == 20) { ans1 -= seg_lnfact(cl); break; }
			else if ((svi = *p) == svim1) ++cl;
			else { rem -= cl; ans1 -= seg_lnfact(cl); if (svi == 0) { ans1 -= seg_lnfact(rem); break; } cl = 1; }
		}
	}
	if (ans1 > -100000.0) { ans2 = seg_lnfact((uint32_t)total); for (int32_t i = 0; sv[i] != 0; ++i) ans2 -= seg_lnfact((uint32_t)sv[i]); }
	return ans1 + ans2 - totseq;
}
void seg_trim(const Letter* seq, int32_t slen, int32_t& leftend, int32_t& rightend) {
	struct Win { int32_t comp[20]{}, sv[21]{}; } base;
	for (int32_t k = 0; k < slen; ++k) if ((uint8_t)seq[k] < 20) seg_inc_sv(base.sv, base.comp[(uint8_t)seq[k]]++);
	int32_t lend = 0, rend = slen - 1, minlen = slen > 51 ? slen - 50 : 1;
	double minprob = 1.0;
	for (int32_t len = slen; len > minlen; --len) {
		Win w = base;
		for (int32_t i = 0;; ++i) {
			const double prob = seg_prob(w.sv, len);
			if (prob < minprob) { minprob = prob; lend = i; rend = len + i - 1; }
			if (i + len >= slen) break;
			if ((uint8_t)seq[i] < 20) seg_dec_sv(w.sv, w.comp[(uint8_t)seq[i]]--);
			if ((uint8_t)seq[i + len] < 20) seg_inc_sv(w.sv, w.comp[(uint8_t)seq[i + len]]++);
		}
		if ((uint8_t)seq[len - 1] < 20) seg_dec_sv(base.sv, base.comp[(uint8_t)seq[len - 1]]--);
	}
	leftend += lend; rightend -= (slen - rend - 1);
}
template<typename Fn> inline void seg_scan(const Letter* seq, int32_t slen, Fn fn) {
	uint8_t comp[21]{}; int S = 0, last = slen - 6;
	for (int32_t j = 0; j < 10; ++j) { const uint8_t c = (uint8_t)seq[j]; if (c < 20) S += kSegDelta[comp[c]++]; else ++comp[20]; }
	for (int32_t i = 4;; ++i) {
		fn(i, comp[20] > 2 || S < kSegHiCut[comp[20]] ? 0 : (S >= kSegLoCut[comp[20]] ? 2 : 1));
		if (i == last) break;
		const uint8_t oc = (uint8_t)seq[i - 4], ic = (uint8_t)seq[i + 6];
		if (oc < 20) S -= kSegDelta[--comp[oc]]; else --comp[20];
		if (ic < 20) S += kSegDelta[comp[ic]++]; else ++comp[20];
	}
}
void seg_sub(const Letter* seq, int32_t slen, int32_t offset, vector<pair<int32_t, int32_t>>& out) {
	if (slen < 10) return;
	uint8_t st_buf[256]; vector<uint8_t> heap_buf;
	uint8_t* st = slen <= 256 ? st_buf : (heap_buf.resize(slen), heap_buf.data());
	seg_scan(seq, slen, [st](int32_t i, uint8_t s) { st[i] = s; });
	for (int32_t i = 4, lowlim = 4, last = slen - 6; i <= last; ++i) if (st[i] == 2) {
		int32_t loi = i, hii = i; while (loi >= lowlim && st[loi] != 0) --loi; while (hii <= last && st[hii] != 0) ++hii;
		int32_t leftend = ++loi - 4, rightend = --hii + 5;
		seg_trim(seq + leftend, rightend - leftend + 1, leftend, rightend);
		if (i + 5 < leftend) {
			const size_t prev = out.size();
			seg_sub(seq + (loi - 4), leftend - (loi - 4), offset + loi - 4, out);
			if (out.size() > prev + 1) { out[prev] = out.back(); out.resize(prev + 1); }
		}
		out.emplace_back(leftend + offset, rightend + offset);
		i = std::min(hii, rightend + 4); lowlim = i + 1;
	}
}
void seg_fast(const Letter* seq, uint32_t len, vector<pair<int32_t, int32_t>>& out) {
	if (len < 10) return;
	int32_t k_l = -1; bool trig = false;
	seg_scan(seq, (int32_t)len, [&](int32_t k, uint8_t s) {
		if (s) { if (k_l < 0) k_l = k; if (s == 2) trig = true; }
		else if (k_l >= 0) { if (trig) seg_sub(seq + (k_l - 4), k - k_l + 9, k_l - 4, out); k_l = -1; trig = false; }
	});
	if (k_l >= 0 && trig) seg_sub(seq + (k_l - 4), (int32_t)len - k_l + 4, k_l - 4, out);
}
}
MaskingStat Masking::operator()(Letter *seq, size_t len, MaskingAlgo algo, const size_t block_id, MaskingTable* table) const
{
	MaskingStat stats;
	if (flag_any(algo, MaskingAlgo::TANTAN)) {
		Mask::Ranges r = Util::tantan::mask(seq, (int)len, (const float**)probMatrixPointersf_, 0.005f, 0.05f, 1.0f / 0.9f, (float)config.tantan_minMaskProb, table ? 0 : 1);
		for (auto i : r) {
			if (table)
				table->add(block_id, i.first, i.second, seq);
			stats.add(MaskingAlgo::TANTAN, i.second - i.first);
		}
	}
	if (flag_any(algo, MaskingAlgo::SEG)) {
		if (blast_seg_->window == 10 && blast_seg_->locut == 1.8 && blast_seg_->hicut == 2.1 && !blast_seg_->overlaps && blast_seg_->maxtrim == 50 && blast_seg_->maxbogus == 2) {
			thread_local vector<pair<int32_t, int32_t>> segs;
			segs.clear(); seg_fast(seq, (uint32_t)len, segs);
			unsigned nMasked = 0;
			for (const auto& s : segs)
				if (table) { table->add(block_id, s.first, s.second + 1, seq); stats.add(MaskingAlgo::SEG, s.second - s.first + 1); }
				else { nMasked += s.second - s.first + 1; std::fill(seq + s.first, seq + s.second + 1, value_traits.mask_char); }
			if (!table) stats.add(MaskingAlgo::SEG, nMasked);
		} else {
		BlastSeqLoc* seg_locs;
		SeqBufferSeg((uint8_t*)seq, (uint32_t)len, 0u, blast_seg_, &seg_locs);
		unsigned nMasked = 0;

		if (seg_locs) {
			BlastSeqLoc* l = seg_locs;
			do {
				if (table) {
					table->add(block_id, l->ssr->left, l->ssr->right + 1, seq);
					stats.add(MaskingAlgo::SEG, l->ssr->right - l->ssr->left + 1);
				}
				else {
					for (signed i = l->ssr->left; i <= l->ssr->right; i++) {
						nMasked++;
						seq[i] = value_traits.mask_char;
					}
				}
			} while ((l = l->next) != 0);
			BlastSeqLocFree(seg_locs);
		}
		if (!table)
			stats.add(MaskingAlgo::SEG, nMasked);
		}
	}
	if (flag_any(algo, MaskingAlgo::MOTIF))
		stats.add(MaskingAlgo::MOTIF, mask_motifs(seq, len, block_id, *table));
	return stats;
}

void Masking::mask_bit(Letter *seq, size_t len) const
{
	Util::tantan::mask(seq, (int)len, (const float**)probMatrixPointersf_, 0.005f, 0.05f, 1.0f / 0.9f, (float)config.tantan_minMaskProb, 2);
}

void Masking::bit_to_hard_mask(Letter *seq, size_t len, size_t &n) const
{
	for (size_t i = 0; i < len; ++i)
		if (seq[i] & bit_mask) {
			seq[i] = value_traits.mask_char;
			++n;
		}
}

void Masking::remove_bit_mask(Letter *seq, size_t len) const
{
	for (size_t i = 0; i < len; ++i)
		if (seq[i] & bit_mask)
			seq[i] &= ~bit_mask;
}

MaskingStat mask_seqs(SequenceSet &seqs, const Masking &masking, bool hard_mask, const MaskingAlgo algo, MaskingTable* table)
{
	MaskingStat stats_all;
	if (algo == MaskingAlgo::NONE)
		return stats_all;
	if (flag_any(algo, MaskingAlgo::MOTIF) && !table)
		throw std::runtime_error("Motif masking requires masking table.");
	vector<thread> threads;
	atomic<BlockId> next(0);	
	std::mutex mtx;

	auto worker = [&]() {
		BlockId i;
		size_t n = 0;
		MaskingStat stats;
		for (BlockId end; (i = next.fetch_add(64, std::memory_order_relaxed)) < seqs.size();)
			for (end = std::min<BlockId>(i + 64, seqs.size()); i < end; ++i)
				if (hard_mask) stats += masking.operator()(seqs.ptr(i), seqs.length(i), algo, i, table); else masking.mask_bit(seqs.ptr(i), seqs.length(i));
		{
			std::lock_guard<std::mutex> lock(mtx);
			stats_all += stats;
		}
		};
	
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(worker);
	for (auto &t : threads)
		t.join();
	return stats_all;
}
