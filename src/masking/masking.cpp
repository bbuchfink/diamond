/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#include <atomic>
#include <numeric>
#include <thread>
#include "masking.h"
#include "tantan/LambdaCalculator.hh"
#include "tantan.h"
#include "blast/blast_filter.h"
#include "data/sequence_set.h"
#include "basic/config.h"

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
	{"seg", MaskingMode::BLAST_SEG}
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
				likelihoodRatioMatrixf_[i][j] = (float)exp(lambda * score_matrix(i, j));
			}
	}
	std::copy(likelihoodRatioMatrixf_, likelihoodRatioMatrixf_ + size, probMatrixPointersf_);

	blast_seg_ = SegParametersNewAa();
}

Masking::~Masking() {
	SegParametersFree(blast_seg_);
}

MaskingStat Masking::operator()(Letter *seq, size_t len, MaskingAlgo algo, const size_t block_id, MaskingTable* table) const
{
	MaskingStat stats;
	if (flag_any(algo, MaskingAlgo::TANTAN)) {
		Mask::Ranges r = Util::tantan::mask(seq, (int)len, (const float**)probMatrixPointersf_, 0.005f, 0.05f, 1.0f / 0.9f, (float)config.tantan_minMaskProb, table ? 0 : 1);
		if (table)
			for (auto i : r) {
				table->add(block_id, i.first, i.second, seq);
				stats.add(MaskingAlgo::TANTAN, i.second - i.first);
			}
	}
	if (flag_any(algo, MaskingAlgo::SEG)) {
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
		while ((i = next.fetch_add(1, std::memory_order_relaxed)) < seqs.size()) {
			if (hard_mask)
				stats += masking.operator()(seqs.ptr(i), seqs.length(i), algo, i, table);
			else
				masking.mask_bit(seqs.ptr(i), seqs.length(i));
		}
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