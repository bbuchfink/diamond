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

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <cstring>
#include <thread>
#include <type_traits>
#include <utility>
#include <numeric>
#include <atomic>
#include <vector>
#include "../dp.h"
#include "swipe.h"
#include "target_iterator.h"
#include "util/data_structures/mem_buffer.h"
#include "util/intrin.h"
#include "util/simd/dispatch.h"
#include "../score_vector_int16.h"
#include "../score_vector_int8.h"

#if !defined(__APPLE__) && !defined(USE_TLS)
#define USE_TLS
#endif

using std::list;
using std::thread;
using std::atomic;
using std::pair;
using std::vector;
using std::runtime_error;

namespace DISPATCH_ARCH {

template<typename Sv>
inline Sv swipe_cell_update(Sv diagonal_cell,
	Sv shift_cell0,
	Sv shift_cell1,
	Sv scores,
	Sv gap_extension,
	Sv gap_open,
	Sv frame_shift,
	Sv &horizontal_gap,
	Sv &vertical_gap)
{
	using std::max;
	const Sv c_diag = diagonal_cell + scores;
	const Sv c_shift = max(shift_cell0, shift_cell1) + (scores - frame_shift);
	const Sv c_gap = max(vertical_gap, horizontal_gap);
	Sv current_cell = max(max(c_diag, c_shift), c_gap);
	saturate(current_cell);
	const Sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap - gap_extension, open);
	horizontal_gap = max(horizontal_gap - gap_extension, open);
	return current_cell;
}

template<typename Sv>
struct Banded3FrameSwipeMatrix
{

	struct ColumnIterator
	{
		ColumnIterator(Sv* hgap_front, Sv* score_front) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front)
		{
			sm4 = ScoreTraits<Sv>::zero();
			sm3 = *score_ptr_;
			sm2 = *(score_ptr_ + 1);
		}
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_;
			sm4 = sm3;
			sm3 = sm2;
			sm2 = *(score_ptr_ + 1);
		}
		inline Sv hgap() const
		{
			return *(hgap_ptr_ + 3);
		}
		inline void set_hgap(const Sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const Sv& x)
		{
			*score_ptr_ = x;
		}
		void set_zero()
		{
			*(score_ptr_ - 1) = ScoreTraits<Sv>::zero();
			*(score_ptr_ - 2) = ScoreTraits<Sv>::zero();
			*(score_ptr_ - 3) = ScoreTraits<Sv>::zero();
		}
		Sv *hgap_ptr_, *score_ptr_;
		Sv sm4, sm3, sm2;
	};

	Banded3FrameSwipeMatrix(size_t band, size_t cols) :
		band_(band)
	{
		hgap_.resize(band + 3);
		score_.resize(band + 1);
		std::fill(hgap_.begin(), hgap_.end(), Sv());
		std::fill(score_.begin(), score_.end(), Sv());
	}

	inline ColumnIterator begin(size_t offset, size_t col)
	{
		return ColumnIterator(&hgap_[offset], &score_[offset]);
	}

	size_t band() const
	{
		return band_;
	}

private:
	const size_t band_;
#ifdef USE_TLS
	static thread_local MemBuffer<Sv> hgap_, score_;
#else
	MemBuffer<Sv> hgap_, score_;
#endif

};

template<typename Sv>
struct Banded3FrameSwipeTracebackMatrix
{

	typedef typename ScoreTraits<Sv>::Score Score;

	struct ColumnIterator
	{
		ColumnIterator(Sv* hgap_front, Sv* score_front, Sv* score_front1) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front),
			score_ptr1_(score_front1)
		{
			sm4 = ScoreTraits<Sv>::zero();
			sm3 = *(score_ptr_++);
			sm2 = *(score_ptr_);
		}
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_; ++score_ptr1_;
			sm4 = sm3;
			sm3 = sm2;
			sm2 = *score_ptr_;
		}
		inline Sv hgap() const
		{
			return *(hgap_ptr_ + 3);
		}
		inline void set_hgap(const Sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const Sv& x)
		{
			*score_ptr1_ = x;
		}
		void set_zero()
		{
			*(score_ptr1_ - 1) = ScoreTraits<Sv>::zero();
			*(score_ptr1_ - 2) = ScoreTraits<Sv>::zero();
			*(score_ptr1_ - 3) = ScoreTraits<Sv>::zero();
		}
		Sv* hgap_ptr_, *score_ptr_, *score_ptr1_;
		Sv sm4, sm3, sm2;
	};

	struct TracebackIterator
	{
		TracebackIterator(const Score *score, size_t band, int frame, int i, int j) :
			band_(band),
			score_(score),
			frame(frame),
			i(i),
			j(j)
		{
			assert(i >= 0 && j >= 0);
		}
		Score score() const
		{
			return *score_;
		}
		Score sm3() const
		{
			return *(score_ - (band_ + 1) * ScoreTraits<Sv>::CHANNELS);
		}
		Score sm4() const
		{
			return *(score_ - (band_ + 2) * ScoreTraits<Sv>::CHANNELS);
		}
		Score sm2() const
		{
			return *(score_ - band_ * ScoreTraits<Sv>::CHANNELS);
		}
		void walk_diagonal()
		{
			score_ -= (band_ + 1) * ScoreTraits<Sv>::CHANNELS;
			--i;
			--j;
			assert(i >= -1 && j >= -1);
		}
		void walk_forward_shift()
		{
			score_ -= (band_ + 2) * ScoreTraits<Sv>::CHANNELS;
			--i;
			--j;
			--frame;
			if (frame == -1) {
				frame = 2;
				--i;
			}
			assert(i >= -1 && j >= -1);
		}
		void walk_reverse_shift()
		{
			score_ -= band_ * ScoreTraits<Sv>::CHANNELS;
			--i;
			--j;
			++frame;
			if (frame == 3) {
				frame = 0;
				++i;
			}
			assert(i >= -1 && j >= -1);
		}
		pair<EditOperation, int> walk_gap(int d0, int d1)
		{
			const int i0 = std::max(d0 + j, 0), j0 = std::max(i - d1, -1);
			const Score *h = score_ - (band_ - 2) * ScoreTraits<Sv>::CHANNELS, *h0 = score_ - (j - j0) * (band_ - 2) * ScoreTraits<Sv>::CHANNELS;
			const Score *v = score_ - 3 * ScoreTraits<Sv>::CHANNELS, *v0 = score_ - (i - i0 + 1) * 3 * ScoreTraits<Sv>::CHANNELS;
			const Score score = this->score();
			const Score e = score_matrix.gap_extend();
			Score g = score_matrix.gap_open() + e;
			int l = 1;
			while (v > v0 && h > h0) {
				if (score + g == *h) {
					walk_hgap(h, l);
					return std::make_pair(op_deletion, l);
				}
				else if (score + g == *v) {
					walk_vgap(v, l);
					return std::make_pair(op_insertion, l);
				}
				h -= (band_ - 2) * ScoreTraits<Sv>::CHANNELS;
				v -= 3 * ScoreTraits<Sv>::CHANNELS;
				++l;
				g += e;
			}
			while (v > v0) {
				if (score + g == *v) {
					walk_vgap(v, l);
					return std::make_pair(op_insertion, l);
				}
				v -= 3 * ScoreTraits<Sv>::CHANNELS;
				++l;
				g += e;
			}
			while (h > h0) {
				if (score + g == *h) {
					walk_hgap(h, l);
					return std::make_pair(op_deletion, l);
				}
				h -= (band_ - 2) * ScoreTraits<Sv>::CHANNELS;
				++l;
				g += e;
			}
			throw runtime_error("Traceback error.");
		}
		void walk_hgap(const Score *h, int l)
		{
			score_ = h;
			j -= l;
			assert(i >= -1 && j >= -1);
		}
		void walk_vgap(const Score *v, int l)
		{
			score_ = v;
			i -= l;
			assert(i >= -1 && j >= -1);
		}
		const size_t band_;
		const Score *score_;
		int frame, i, j;
	};

	TracebackIterator traceback(size_t col, int i0, int j, int dna_len, size_t channel, Score score) const
	{
		const int i_ = std::max(-i0, 0) * 3,
			i1 = (int)std::min(band_, size_t(dna_len - 2 - i0 * 3));
		const Score *s = (Score*)(&score_[col*(band_ + 1) + i_]) + channel;
		for (int i = i_; i < i1; ++i, s += ScoreTraits<Sv>::CHANNELS)
			if (*s == score)
				return TracebackIterator(s, band_, i % 3, i0 + i / 3, j);
		throw runtime_error("Trackback error.");
	}

	Banded3FrameSwipeTracebackMatrix(size_t band, size_t cols) :
		band_(band)
	{
		hgap_.resize(band + 3);
		score_.resize((band + 1) * (cols + 1));
		const Sv z = Sv();
		std::fill(hgap_.begin(), hgap_.end(), z);
		std::fill(score_.begin(), score_.begin() + band + 1, z);
		for (size_t i = 0; i < cols; ++i)
			score_[i*(band + 1) + band] = z;
	}

	inline ColumnIterator begin(size_t offset, size_t col)
	{
		return ColumnIterator(&hgap_[offset], &score_[col*(band_ + 1) + offset], &score_[(col + 1)*(band_ + 1) + offset]);
	}

	size_t band() const
	{
		return band_;
	}

private:
	const size_t band_;
#ifdef USE_TLS
	static thread_local MemBuffer<Sv> hgap_, score_;
#else
	MemBuffer<Sv> hgap_, score_;
#endif

};

#ifdef USE_TLS
template<typename Sv> thread_local MemBuffer<Sv> Banded3FrameSwipeMatrix<Sv>::hgap_;
template<typename Sv> thread_local MemBuffer<Sv> Banded3FrameSwipeMatrix<Sv>::score_;
template<typename Sv> thread_local MemBuffer<Sv> Banded3FrameSwipeTracebackMatrix<Sv>::hgap_;
template<typename Sv> thread_local MemBuffer<Sv> Banded3FrameSwipeTracebackMatrix<Sv>::score_;
#endif

template<typename Sv, typename Traceback>
struct Banded3FrameSwipeMatrixRef
{
};

template<typename Sv>
struct Banded3FrameSwipeMatrixRef<Sv, DP::Traceback>
{
	typedef Banded3FrameSwipeTracebackMatrix<Sv> type;
};

template<typename Sv>
struct Banded3FrameSwipeMatrixRef<Sv, DP::ScoreOnly>
{
	typedef Banded3FrameSwipeMatrix<Sv> type;
};

template<typename Sv>
Hsp traceback(Sequence *query, Strand strand, int dna_len, const Banded3FrameSwipeTracebackMatrix<Sv> &dp, const DpTarget &target, int d_begin,
	typename ScoreTraits<Sv	>::Score max_score, double evalue, int max_col, int channel, int i0, int i1)
{
	typedef typename ScoreTraits<Sv>::Score Score;
	const int j0 = i1 - (target.d_end - 1), d1 = target.d_end;
	typename Banded3FrameSwipeTracebackMatrix<Sv>::TracebackIterator it(dp.traceback(max_col + 1, i0 + max_col, j0 + max_col, dna_len, channel, max_score));
	
	Hsp out(true);
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<Sv>::int_score(max_score) * config.cbs_matrix_scale;
	out.bit_score = score_matrix.bitscore(out.score);
	out.evalue = evalue;
	out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));

	out.set_end(it.i + 1, it.j + 1, Frame(strand, it.frame), dna_len);

	while (it.score() > ScoreTraits<Sv>::zero_score()) {
		const Letter q = query[it.frame][it.i], s = target.seq[it.j];
		const Score m = score_matrix(q, s), score = it.score();
		if (score == it.sm3() + m) {
			out.push_match(q, s, m > (Score)0);
			it.walk_diagonal();
		}
		else if (score == it.sm4() + m - score_matrix.frame_shift()) {
			out.push_match(q, s, m > (Score)0);
			out.transcript.push_back(op_frameshift_forward);
			it.walk_forward_shift();
		}
		else if (score == it.sm2() + m - score_matrix.frame_shift()) {
			out.push_match(q, s, m > (Score)0);
			out.transcript.push_back(op_frameshift_reverse);
			it.walk_reverse_shift();
		}
		else {
			const pair<EditOperation, int> g(it.walk_gap(d_begin, d1));
			out.push_gap(g.first, g.second, target.seq.data() + it.j + g.second);
		}
	}

	out.set_begin(it.i + 1, it.j + 1, Frame(strand, it.frame), dna_len);
	out.transcript.reverse();
	out.transcript.push_terminator();
	out.subject_source_range = out.subject_range;
	return out;
}

template<typename Sv>
Hsp traceback(Sequence *query, Strand strand, int dna_len, const Banded3FrameSwipeMatrix<Sv> &dp, const DpTarget &target, int d_begin,
	typename ScoreTraits<Sv>::Score max_score, double evalue, int max_col, int channel, int i0, int i1)
{
	Hsp out(false);
	const int j0 = i1 - (target.d_end - 1);
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<Sv>::int_score(max_score) * config.cbs_matrix_scale;
	out.bit_score = score_matrix.bitscore(out.score);
	out.evalue = evalue;
	out.query_range.end_ = std::min(i0 + max_col + (int)dp.band() / 3 / 2, (int)query[0].length());
	out.query_range.begin_ = std::max(out.query_range.end_ - (j0 + max_col), 0);
	out.frame = strand == FORWARD ? 0 : 3;
	out.query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(out.query_range.begin_, Frame(out.frame)), TranslatedPosition(out.query_range.end_, Frame(out.frame)), dna_len);
	return out;
}

template<typename Sv, typename Traceback>
list<Hsp> banded_3frame_swipe(
	const TranslatedSequence& query,
	Strand strand, std::vector<DpTarget>::const_iterator subject_begin,
	std::vector<DpTarget>::const_iterator subject_end,
	DpStat &stat,
	bool parallel,
	std::vector<DpTarget> &overflow)
{
	typedef typename Banded3FrameSwipeMatrixRef<Sv, Traceback>::type Matrix;
	typedef typename ScoreTraits<Sv>::Score Score;

	assert(subject_end - subject_begin <= ScoreTraits<Sv>::CHANNELS);
	Sequence q[3];
	query.get_strand(strand, q);
	const int qlen = (int)q[0].length(), qlen2 = (int)q[1].length(), qlen3 = (int)q[2].length();

	int band = 0;
	for (vector<DpTarget>::const_iterator j = subject_begin; j < subject_end; ++j)
		band = std::max(band, j->d_end - j->d_begin);

	using std::max;
	int i0 = INT_MAX, i1 = INT_MAX, d_begin[ScoreTraits<Sv>::CHANNELS];
	for (vector<DpTarget>::const_iterator j = subject_begin; j < subject_end; ++j) {
		d_begin[j - subject_begin] = j->d_end - band;
		int i2 = std::max(j->d_end - 1, 0);
		i1 = std::min(i1, i2);
		i0 = std::min(i0, i2 + 1 - band);
	}
	const int orig_i0 = i0, orig_i1 = i1;

	TargetIterator<Sv> targets(subject_begin, subject_end, false, i1, qlen, d_begin);
	Matrix dp(band * 3, targets.cols);

	const Sv open_penalty(score_matrix.gap_open() + score_matrix.gap_extend()),
		extend_penalty(score_matrix.gap_extend()),
		frameshift_penalty(score_matrix.frame_shift());
	
	SwipeProfile<Sv> profile;
#ifndef __SSSE3__
	std::array<const int8_t*, 32> target_scores;
#endif
	Score best[ScoreTraits<Sv>::CHANNELS];
	int max_col[ScoreTraits<Sv>::CHANNELS];
	for (int i = 0; i < ScoreTraits<Sv>::CHANNELS; ++i) {
		best[i] = ScoreTraits<Sv>::zero_score();
		max_col[i] = 0;
	}

	const int xdrop = config.raw_swipe_xdrop;

	int j = 0;
	while (targets.active.size() > 0) {
		const int i0_ = std::max(i0, 0), i1_ = std::min(i1, qlen - 1);
		if (i0_ > i1_)
			break;
		typename Matrix::ColumnIterator it(dp.begin((i0_ - i0) * 3, j));
		if (i0_ - i0 > 0)
			it.set_zero();
		Sv vgap0, vgap1, vgap2, hgap0, hgap1, hgap2;
		vgap0 = vgap1 = vgap2 = ScoreTraits<Sv>::zero();
		Sv col_best0, col_best1, col_best2, col_best;
		col_best0 = col_best1 = col_best2 = col_best = ScoreTraits<Sv>::zero();

#ifdef __SSSE3__
		profile.set(targets.get());
#else
		profile.set(targets.get(target_scores.data()));
#endif
		const int common_i1 = std::min(i1_, std::min(qlen2, qlen3) - 1);
		int i = i0_;
		for (; i <= common_i1; ++i) {
			const Sv s0 = profile.get(letter_mask(q[0][i]));
			const Sv s1 = profile.get(letter_mask(q[1][i]));
			const Sv s2 = profile.get(letter_mask(q[2][i]));

			hgap0 = it.hgap();
			Sv next0 = swipe_cell_update<Sv>(it.sm3, it.sm4, it.sm2, s0, extend_penalty, open_penalty, frameshift_penalty, hgap0, vgap0);
			it.set_hgap(hgap0);
			it.set_score(next0);
			col_best0 = max(col_best0, next0);
			++it;

			hgap1 = it.hgap();
			Sv next1 = swipe_cell_update<Sv>(it.sm3, it.sm4, it.sm2, s1, extend_penalty, open_penalty, frameshift_penalty, hgap1, vgap1);
			it.set_hgap(hgap1);
			it.set_score(next1);
			col_best1 = max(col_best1, next1);
			++it;

			hgap2 = it.hgap();
			Sv next2 = swipe_cell_update<Sv>(it.sm3, it.sm4, it.sm2, s2, extend_penalty, open_penalty, frameshift_penalty, hgap2, vgap2);
			it.set_hgap(hgap2);
			it.set_score(next2);
			col_best2 = max(col_best2, next2);
			++it;
		}

		col_best = max(col_best0, max(col_best1, col_best2));

		for (; i <= i1_; ++i) {
			hgap0 = it.hgap();
			Sv next0 = swipe_cell_update<Sv>(it.sm3, it.sm4, it.sm2, profile.get(letter_mask(q[0][i])), extend_penalty, open_penalty, frameshift_penalty, hgap0, vgap0);
			it.set_hgap(hgap0);
			it.set_score(next0);
			col_best = max(col_best, next0);
			++it;

			if (i >= qlen2)
				break;
			hgap1 = it.hgap();
			Sv next1 = swipe_cell_update<Sv>(it.sm3, it.sm4, it.sm2, profile.get(letter_mask(q[1][i])), extend_penalty, open_penalty, frameshift_penalty, hgap1, vgap1);
			it.set_hgap(hgap1);
			it.set_score(next1);
			col_best = max(col_best, next1);
			++it;

			if (i >= qlen3)
				break;
			hgap2 = it.hgap();
			Sv next2 = swipe_cell_update<Sv>(it.sm3, it.sm4, it.sm2, profile.get(letter_mask(q[2][i])), extend_penalty, open_penalty, frameshift_penalty, hgap2, vgap2);
			it.set_hgap(hgap2);
			it.set_score(next2);
			col_best = max(col_best, next2);
			++it;
		}

#ifdef DP_STAT
		//stat.net_cells += targets.live * (i1_ - i0_ + 1) * 3;
		//stat.gross_cells += ScoreTraits<_sv>::CHANNELS * (i1_ - i0_ + 1) * 3;
#endif

		Score col_best_[ScoreTraits<Sv>::CHANNELS];
		store_sv(col_best, col_best_);
		if (xdrop > 0) {
			for (int i = 0; i < targets.active.size();) {
				const int channel = targets.active[i];
				if (col_best_[channel] > best[channel]) {
					best[channel] = col_best_[channel];
					max_col[channel] = j;
				}
				if (best[channel] - col_best_[channel] >= xdrop || !targets.inc(channel))
					targets.active.erase(i);
				else
					++i;
			}
		}
		else {
			for (int i = 0; i < targets.active.size();) {
				int channel = targets.active[i];
				if (!targets.inc(channel))
					targets.active.erase(i);
				else
					++i;
				if (col_best_[channel] > best[channel]) {
					best[channel] = col_best_[channel];
					max_col[channel] = j;
				}
			}
		}
		++i0;
		++i1;
		++j;
	}
	
	list<Hsp> out;
	for (int i = 0; i < targets.n_targets; ++i) {
		if (best[i] < ScoreTraits<Sv>::max_score()) {
			const int score = ScoreTraits<Sv>::int_score(best[i]) * config.cbs_matrix_scale;
			const double evalue = score_matrix.evalue(score, qlen, (unsigned)subject_begin[i].seq.length());
			if (score_matrix.report_cutoff(score, evalue))
				out.push_back(traceback<Sv>(q, strand, (int)query.source().length(), dp, subject_begin[i], d_begin[i], best[i], evalue, max_col[i], i, orig_i0, orig_i1));
		}
		else
			overflow.push_back(subject_begin[i]);
	}
	return out;
}

template<typename Sv>
list<Hsp> banded_3frame_swipe_targets(std::vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	bool score_only,
	const TranslatedSequence &query,
	Strand strand,
	DpStat &stat,
	bool parallel,
	std::vector<DpTarget> &overflow)
{
	list<Hsp> out;
	for (vector<DpTarget>::const_iterator i = begin; i < end; i += std::min((ptrdiff_t)ScoreTraits<Sv>::CHANNELS, end - i)) {
		if (score_only)
			out.splice(out.end(), banded_3frame_swipe<Sv, DP::ScoreOnly>(query, strand, i, i + std::min(ptrdiff_t(ScoreTraits<Sv>::CHANNELS), end - i), stat, parallel, overflow));
		else
			out.splice(out.end(), banded_3frame_swipe<Sv, DP::Traceback>(query, strand, i, i + std::min(ptrdiff_t(ScoreTraits<Sv>::CHANNELS), end - i), stat, parallel, overflow));
	}
	return out;
}

void banded_3frame_swipe_worker(std::vector<DpTarget>::const_iterator begin,
	std::vector<DpTarget>::const_iterator end,
	atomic<size_t> *next,
	bool score_only,
	const TranslatedSequence *query,
	Strand strand,
	list<Hsp> *out,
	vector<DpTarget> *overflow)
{
	DpStat stat;
	size_t pos;
	vector<DpTarget> of;
	while (begin + (pos = next->fetch_add(config.swipe_chunk_size, std::memory_order_relaxed)) < end)
#if defined(__SSE2__) | defined(__ARM_NEON)
		if(score_only)
			out->splice(out->end(), banded_3frame_swipe_targets<ScoreVector<int16_t, SHRT_MIN>>(begin + pos, min(begin + pos + config.swipe_chunk_size, end), score_only, *query, strand, stat, true, of));
		else
			out->splice(out->end(), banded_3frame_swipe_targets<int32_t>(begin + pos, min(begin + pos + config.swipe_chunk_size, end), score_only, *query, strand, stat, true, of));
#else
		out->splice(out->end(), banded_3frame_swipe_targets<int32_t>(begin + pos, min(begin + pos + config.swipe_chunk_size, end), score_only, *query, strand, stat, true, of));
#endif
	*overflow = std::move(of);
}

#if defined(__AVX2__) && ARCH_AVX2_KERNELS
static constexpr int ANTIDIAG_SWIPE_WEIGHT = 4;
struct alignas(16) AntidiagSeam { int16_t s0, s1, s2, h0, h1, h2, p0, p1; };
struct AntidiagStripMeta { int i_base, s_end; size_t dir_off; };
struct AntidiagCluster { int q0, q1, pstride; size_t prof_off; };
struct AntidiagPass1Hit { const Letter* seq; int32_t d_begin, d_end, best_i, best_j; int16_t best_sc; };
struct AntidiagCache {
	const int8_t* mat = nullptr;
	vector<Letter> q_copy[3];
	vector<AntidiagCluster> clusters;
	vector<int8_t> prof_buf;
	size_t prof_used = 0;
	vector<AntidiagPass1Hit> hits;
};
struct AntidiagWorkspace {
	vector<AntidiagSeam> seam;
	vector<uint8_t> dir;
	vector<AntidiagStripMeta> strips;
	AntidiagCache cache[2];
	void trim(bool score_only) {
		if (score_only) return;
		if (seam.capacity() > (1u << 18)) vector<AntidiagSeam>().swap(seam);
		if (dir.capacity() > (1u << 24)) vector<uint8_t>().swap(dir);
		for (int s = 0; s < 2; ++s) {
			if (cache[s].prof_buf.capacity() > (1u << 24)) {
				vector<int8_t>().swap(cache[s].prof_buf);
				cache[s].clusters.clear();
				cache[s].prof_used = 0;
				cache[s].hits.clear();
				cache[s].mat = nullptr;
			}
		}
	}
};

FORCE_INLINE __m256i antidiag_shift_zero(__m256i v) { return _mm256_alignr_epi8(v, _mm256_permute2x128_si256(v, v, 0x08), 14); }
FORCE_INLINE __m256i antidiag_shift_lane0(__m256i v, int16_t x) { return _mm256_alignr_epi8(v, _mm256_permute2x128_si256(_mm256_set1_epi16(x), v, 0x20), 14); }
FORCE_INLINE void antidiag_store_seam(AntidiagSeam* dst, __m256i s0, __m256i s1, __m256i s2, __m256i h0, __m256i h1, __m256i h2) {
	const __m256i p01 = _mm256_unpackhi_epi16(s0, s1), p23 = _mm256_unpackhi_epi16(s2, h0), p45 = _mm256_unpackhi_epi16(h1, h2);
	const __m256i r = _mm256_unpackhi_epi64(_mm256_unpackhi_epi32(p01, p23), _mm256_unpackhi_epi32(p45, p45));
	_mm_storeu_si128((__m128i*)dst, _mm256_extracti128_si256(r, 1));
}
FORCE_INLINE __m256i antidiag_seed_lane0(int16_t lane0) {
	return _mm256_inserti128_si256(_mm256_setzero_si256(), _mm_cvtsi32_si128((int)(uint16_t)lane0), 0);
}
FORCE_INLINE __m256i antidiag_lane_mask(int i0, int qlen_f, int strip_w) {
	const __m256i idx = _mm256_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
	const int hi = std::min(strip_w, std::min(i0 + 1, 16)), lo = std::max(-1, std::min(i0 - qlen_f, 16));
	return _mm256_and_si256(_mm256_cmpgt_epi16(_mm256_set1_epi16((int16_t)hi), idx),
		_mm256_cmpgt_epi16(idx, _mm256_set1_epi16((int16_t)lo)));
}

FORCE_INLINE void antidiag_transpose_f01(const signed char* const* base_p, int b, int fstride, int sw, int blk_len, __m256i* dst) {
	__m256i r[16];
	if (sw == 16) {
		for (int u = 0; u < 16; ++u) {
			const signed char* p = base_p[u] + b;
			r[u] = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i*)p)), _mm_loadu_si128((const __m128i*)(p + fstride)), 1);
		}
	} else {
		for (int u = 0; u < sw; ++u) {
			const signed char* p = base_p[u] + b;
			r[u] = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i*)p)), _mm_loadu_si128((const __m128i*)(p + fstride)), 1);
		}
		for (int u = sw; u < 16; ++u) r[u] = _mm256_setzero_si256();
	}
	__m256i t[16], u[16], v[16];
	for (int i = 0; i < 8; ++i) {
		t[2 * i] = _mm256_unpacklo_epi8(r[2 * i], r[2 * i + 1]);
		t[2 * i + 1] = _mm256_unpackhi_epi8(r[2 * i], r[2 * i + 1]);
	}
	for (int g = 0; g < 4; ++g) {
		const int o = g * 4;
		u[o] = _mm256_unpacklo_epi16(t[o], t[o + 2]); u[o + 1] = _mm256_unpackhi_epi16(t[o], t[o + 2]);
		u[o + 2] = _mm256_unpacklo_epi16(t[o + 1], t[o + 3]); u[o + 3] = _mm256_unpackhi_epi16(t[o + 1], t[o + 3]);
	}
	for (int i = 0; i < 4; ++i) {
		v[2 * i] = _mm256_unpacklo_epi32(u[i], u[i + 4]); v[2 * i + 1] = _mm256_unpackhi_epi32(u[i], u[i + 4]);
		v[2 * i + 8] = _mm256_unpacklo_epi32(u[i + 8], u[i + 12]); v[2 * i + 9] = _mm256_unpackhi_epi32(u[i + 8], u[i + 12]);
	}
	const int pairs = (blk_len <= 8) ? 4 : 8;
	for (int i = 0; i < pairs; ++i) {
		_mm256_store_si256(&dst[2 * i], _mm256_unpacklo_epi64(v[i], v[i + 8]));
		_mm256_store_si256(&dst[2 * i + 1], _mm256_unpackhi_epi64(v[i], v[i + 8]));
	}
}

FORCE_INLINE void antidiag_transpose_f2(const signed char* const* base_p, int b, int fstride2, int sw, int blk_len, __m128i* dst) {
	__m128i r[16];
	if (sw == 16) {
		for (int u = 0; u < 16; ++u) r[u] = _mm_loadu_si128((const __m128i*)(base_p[u] + fstride2 + b));
	} else {
		for (int u = 0; u < sw; ++u) r[u] = _mm_loadu_si128((const __m128i*)(base_p[u] + fstride2 + b));
		for (int u = sw; u < 16; ++u) r[u] = _mm_setzero_si128();
	}
	__m128i t[16], u[16], v[16];
	for (int i = 0; i < 8; ++i) {
		t[2 * i] = _mm_unpacklo_epi8(r[2 * i], r[2 * i + 1]);
		t[2 * i + 1] = _mm_unpackhi_epi8(r[2 * i], r[2 * i + 1]);
	}
	for (int g = 0; g < 4; ++g) {
		const int o = g * 4;
		u[o] = _mm_unpacklo_epi16(t[o], t[o + 2]); u[o + 1] = _mm_unpackhi_epi16(t[o], t[o + 2]);
		u[o + 2] = _mm_unpacklo_epi16(t[o + 1], t[o + 3]); u[o + 3] = _mm_unpackhi_epi16(t[o + 1], t[o + 3]);
	}
	for (int i = 0; i < 4; ++i) {
		v[2 * i] = _mm_unpacklo_epi32(u[i], u[i + 4]); v[2 * i + 1] = _mm_unpackhi_epi32(u[i], u[i + 4]);
		v[2 * i + 8] = _mm_unpacklo_epi32(u[i + 8], u[i + 12]); v[2 * i + 9] = _mm_unpackhi_epi32(u[i + 8], u[i + 12]);
	}
	const int pairs = (blk_len <= 8) ? 4 : 8;
	for (int i = 0; i < pairs; ++i) {
		_mm_store_si128(&dst[2 * i], _mm_unpacklo_epi64(v[i], v[i + 8]));
		_mm_store_si128(&dst[2 * i + 1], _mm_unpackhi_epi64(v[i], v[i + 8]));
	}
}

template<bool Traceback, bool Masked>
FORCE_INLINE __m256i antidiag_step_frame(__m256i sm3, __m256i sm4, __m256i sm2, __m256i m,
	__m256i v_fs, __m256i v_open, __m256i v_ext, __m256i mask, __m256i& v, __m256i& h, __m256i& dir_out)
{
	const __m256i sm42 = _mm256_max_epi16(sm4, sm2), fs = _mm256_subs_epi16(sm42, v_fs);
	const __m256i diag = _mm256_adds_epi16(_mm256_max_epi16(sm3, fs), m);
	const __m256i h_in = h, vh = _mm256_max_epi16(v, h_in);
	__m256i cur = _mm256_max_epi16(diag, vh);
	if (Masked) cur = _mm256_and_si256(cur, mask);
	const __m256i open = _mm256_subs_epu16(cur, v_open), ve = _mm256_subs_epu16(v, v_ext), he = _mm256_subs_epu16(h_in, v_ext);
	if (Traceback) {
		const __m256i is_vh = _mm256_cmpgt_epi16(vh, diag);
		const __m256i b2 = _mm256_blendv_epi8(_mm256_cmpgt_epi16(sm2, sm4), _mm256_cmpgt_epi16(v, h_in), is_vh);
		__m256i lo = _mm256_srli_epi16(b2, 15), hi = _mm256_srli_epi16(_mm256_cmpgt_epi16(ve, open), 15);
		lo = _mm256_sub_epi16(_mm256_add_epi16(lo, lo), _mm256_cmpgt_epi16(fs, sm3));
		lo = _mm256_sub_epi16(_mm256_add_epi16(lo, lo), is_vh);
		hi = _mm256_sub_epi16(_mm256_add_epi16(hi, hi), _mm256_cmpgt_epi16(he, open));
		dir_out = _mm256_or_si256(_mm256_slli_epi16(hi, 3), lo);
	}
	v = _mm256_max_epi16(ve, open);
	h = _mm256_max_epi16(he, open);
	if (Masked) { v = _mm256_and_si256(v, mask); h = _mm256_and_si256(h, mask); }
	return cur;
}

template<bool Traceback>
static int antidiag_run_target_avx2(const Sequence* q, const int8_t* prof, int pstride, int q0,
	const DpTarget& t, int64_t j_beg, int64_t j_end, int i_max, int16_t init_best, int xdrop,
	AntidiagWorkspace& ws, int& best_i, int& best_j, int& best_f)
{
	const int full_qlen = (int)q[0].length();
	const int qlen = std::min(full_qlen, i_max + 1), qlen2 = std::min((int)q[1].length(), i_max + 1), qlen3 = std::min((int)q[2].length(), i_max + 1);
	AntidiagSeam* seam = nullptr;
	if (j_end - j_beg > 16) {
		if ((int)ws.seam.size() < full_qlen + 48) ws.seam.resize((size_t)full_qlen + 48);
		seam = ws.seam.data() + 16;
	}
	int prev_w_lo = 0, prev_w_hi = -1;
	if (Traceback) {
		size_t total_dir = 0;
		ws.strips.clear();
		for (int64_t j0 = j_beg; j0 < j_end; j0 += 16) {
			const int sw = (int)std::min<int64_t>(16, j_end - j0);
			const int ib = std::max(0, (int)j0 + t.d_begin);
			const int se = std::max(0, std::min(qlen - 1, (int)(j0 + sw - 1) + t.d_end - 1) - ib + sw);
			ws.strips.push_back({ ib, se, total_dir });
			total_dir += (size_t)se * 32;
		}
		if (ws.dir.size() < total_dir) ws.dir.resize(total_dir);
	}
	const __m256i v_open = _mm256_set1_epi16((int16_t)(score_matrix.gap_open() + score_matrix.gap_extend()));
	const __m256i v_ext = _mm256_set1_epi16((int16_t)score_matrix.gap_extend());
	const __m256i v_fs = _mm256_set1_epi16((int16_t)score_matrix.frame_shift());
	const int fstride = pstride, fstride2 = 2 * pstride;
	int16_t best_sc = init_best;
	best_i = -1;
	int xdrop_best = 0;
	__m256i v_best = _mm256_set1_epi16(best_sc);
	alignas(32) __m256i trans01[16];
	alignas(16) __m128i trans2[16];
	const Letter* tseq = t.seq.data();
	size_t strip_idx = 0;
	for (int64_t j0 = j_beg; j0 < j_end; j0 += 16, ++strip_idx) {
		const int sw = (int)std::min<int64_t>(16, j_end - j0);
		const int i_base = std::max(0, (int)j0 + t.d_begin);
		const int s_end = std::min(qlen - 1, (int)(j0 + sw - 1) + t.d_end - 1) - i_base + sw;
		if (s_end <= 0) {
			prev_w_lo = 0; prev_w_hi = -1;
			if (xdrop > 0 && xdrop_best >= xdrop) break;
			continue;
		}
		if (j0 > j_beg) {
			const int r_lo = i_base - 2, r_hi = i_base + s_end;
			if (prev_w_hi < prev_w_lo) {
				std::fill(seam + r_lo, seam + r_hi + 1, AntidiagSeam{});
			} else {
				if (r_lo < prev_w_lo) std::fill(seam + r_lo, seam + std::min(r_hi, prev_w_lo - 1) + 1, AntidiagSeam{});
				if (r_hi > prev_w_hi) std::fill(seam + std::max(r_lo, prev_w_hi + 1), seam + r_hi + 1, AntidiagSeam{});
			}
		}
		uint8_t* dir_base = Traceback ? (ws.dir.data() + ws.strips[strip_idx].dir_off) : nullptr;
		const signed char* base_p[16];
		for (int u = 0; u < sw; ++u) {
			const unsigned c = (unsigned)letter_mask(tseq[j0 + u]);
			const int off = (int)(c < AMINO_ACID_COUNT ? c : (unsigned)MASK_LETTER) * 3 * pstride + 32 + (i_base - u - q0);
			base_p[u] = (const signed char*)prof + off;
		}
		__m256i strip_max_v = _mm256_setzero_si256();
		auto run_strip = [&](auto hl_tag, auto hr_tag) {
			constexpr bool HasL = decltype(hl_tag)::value, HasR = decltype(hr_tag)::value;
			__m256i v0 = _mm256_setzero_si256(), v1 = v0, v2 = v0;
			__m256i h0 = HasL ? antidiag_seed_lane0(seam[i_base].h0) : v0, h1 = HasL ? antidiag_seed_lane0(seam[i_base].h1) : v0, h2 = HasL ? antidiag_seed_lane0(seam[i_base].h2) : v0;
			__m256i d1_0 = HasL ? antidiag_seed_lane0(seam[i_base].s0) : v0, d1_1 = HasL ? antidiag_seed_lane0(seam[i_base].s1) : v0, d1_2 = HasL ? antidiag_seed_lane0(seam[i_base].s2) : v0;
			__m256i d2_0 = HasL ? antidiag_seed_lane0(seam[i_base - 1].s0) : v0, d2_1 = HasL ? antidiag_seed_lane0(seam[i_base - 1].s1) : v0, d2_2 = HasL ? antidiag_seed_lane0(seam[i_base - 1].s2) : v0;
			__m256i d3_2 = HasL ? antidiag_seed_lane0(seam[i_base - 2].s2) : v0;
			for (int b = 0; b < s_end; b += 16) {
				const int i_blk = i_base + b, blk_len = std::min(16, s_end - b);
				antidiag_transpose_f01(base_p, b, fstride, sw, blk_len, trans01);
				antidiag_transpose_f2(base_p, b, fstride2, sw, blk_len, trans2);
				auto step_k = [&](auto m_tag, int k)
#if defined(__GNUC__)
					__attribute__((always_inline))
#endif
				{
					constexpr bool M = decltype(m_tag)::value;
					const int i0 = i_blk + k;
					__m256i dir0, dir1, dir2;
					const __m256i mk0 = M ? antidiag_lane_mask(i0, qlen, sw) : _mm256_setzero_si256();
					const __m256i m0 = _mm256_cvtepi8_epi16(_mm_load_si128((const __m128i*)&trans01[k]));
					const __m256i s0 = antidiag_step_frame<Traceback, M>(d2_0, d3_2, d2_1, m0, v_fs, v_open, v_ext, mk0, v0, h0, dir0);
					const __m256i mk1 = (!M || i0 < qlen3 || qlen2 == qlen) ? mk0 : antidiag_lane_mask(i0, qlen2, sw);
					const __m256i m1 = _mm256_cvtepi8_epi16(_mm_load_si128((const __m128i*)&trans01[k] + 1));
					const __m256i s1 = antidiag_step_frame<Traceback, M>(d2_1, d2_0, d2_2, m1, v_fs, v_open, v_ext, mk1, v1, h1, dir1);
					d3_2 = d2_2; d2_0 = d1_0;
					const __m256i mk2 = (!M || i0 < qlen3 || qlen3 == qlen2) ? mk1 : antidiag_lane_mask(i0, qlen3, sw);
					const __m256i m2 = _mm256_cvtepi8_epi16(_mm_load_si128(&trans2[k]));
					const __m256i s2 = antidiag_step_frame<Traceback, M>(d3_2, d2_1, d2_0, m2, v_fs, v_open, v_ext, mk2, v2, h2, dir2);
					if (Traceback) {
						const __m256i dir3 = _mm256_or_si256(dir0, _mm256_or_si256(_mm256_slli_epi16(dir1, 5), _mm256_slli_epi16(dir2, 10)));
						_mm256_storeu_si256((__m256i*)(dir_base + (size_t)(b + k) * 32), dir3);
					}
					d2_1 = d1_1; d2_2 = d1_2;
					const __m256i smax = _mm256_max_epi16(_mm256_max_epi16(s0, s1), s2);
					if (xdrop > 0) strip_max_v = _mm256_max_epi16(strip_max_v, smax);
					const int msk = _mm256_movemask_epi8(_mm256_cmpgt_epi16(smax, v_best));
					if (UNLIKELY(msk != 0)) {
						if (!Traceback) {
							alignas(32) int16_t bm[16]; _mm256_store_si256((__m256i*)bm, smax);
							for (uint32_t m = (uint32_t)msk; m; ) {
								const int u = ctz(m) >> 1; m &= ~(3u << (u * 2));
								if (bm[u] > best_sc) { best_sc = bm[u]; best_i = i0 - u; best_j = (int)j0 + u; }
							}
						} else {
							alignas(32) int16_t b0[16], b1[16], b2[16];
							_mm256_store_si256((__m256i*)b0, s0); _mm256_store_si256((__m256i*)b1, s1); _mm256_store_si256((__m256i*)b2, s2);
							for (uint32_t m = (uint32_t)msk; m; ) {
								const int u = ctz(m) >> 1; m &= ~(3u << (u * 2));
								if (b0[u] > best_sc) { best_sc = b0[u]; best_i = i0 - u; best_j = (int)j0 + u; best_f = 0; }
								if (b1[u] > best_sc) { best_sc = b1[u]; best_i = i0 - u; best_j = (int)j0 + u; best_f = 1; }
								if (b2[u] > best_sc) { best_sc = b2[u]; best_i = i0 - u; best_j = (int)j0 + u; best_f = 2; }
							}
						}
						v_best = _mm256_set1_epi16(best_sc);
					}
					if (HasR) antidiag_store_seam(&seam[i0 - 15], s0, s1, s2, h0, h1, h2);
					if (HasL) {
						const AntidiagSeam& nx = seam[i0 + 1];
						d1_0 = antidiag_shift_lane0(s0, nx.s0); d1_1 = antidiag_shift_lane0(s1, nx.s1); d1_2 = antidiag_shift_lane0(s2, nx.s2);
						h0 = antidiag_shift_lane0(h0, nx.h0); h1 = antidiag_shift_lane0(h1, nx.h1); h2 = antidiag_shift_lane0(h2, nx.h2);
					} else {
						d1_0 = antidiag_shift_zero(s0); d1_1 = antidiag_shift_zero(s1); d1_2 = antidiag_shift_zero(s2);
						h0 = antidiag_shift_zero(h0); h1 = antidiag_shift_zero(h1); h2 = antidiag_shift_zero(h2);
					}
				};
				if (sw == 16 && blk_len == 16 && i_blk >= 15 && i_blk + 15 < qlen3)
					for (int k = 0; k < 16; ++k) step_k(std::false_type{}, k);
				else
					for (int k = 0; k < blk_len; ++k) step_k(std::true_type{}, k);
			}
		};
		if (j0 > j_beg) { if (j0 + 16 < j_end) run_strip(std::true_type{}, std::true_type{}); else run_strip(std::true_type{}, std::false_type{}); }
		else { if (j0 + 16 < j_end) run_strip(std::false_type{}, std::true_type{}); else run_strip(std::false_type{}, std::false_type{}); }
		prev_w_lo = i_base - 15; prev_w_hi = i_base + s_end - 16;
		if (xdrop > 0) {
			__m128i sm128 = _mm_max_epi16(_mm256_castsi256_si128(strip_max_v), _mm256_extracti128_si256(strip_max_v, 1));
			sm128 = _mm_max_epi16(sm128, _mm_shuffle_epi32(sm128, _MM_SHUFFLE(1, 0, 3, 2)));
			sm128 = _mm_max_epi16(sm128, _mm_shuffle_epi32(sm128, _MM_SHUFFLE(2, 3, 0, 1)));
			sm128 = _mm_max_epi16(sm128, _mm_shufflelo_epi16(sm128, _MM_SHUFFLE(2, 3, 0, 1)));
			const int strip_max = (int)(int16_t)_mm_extract_epi16(sm128, 0);
			if (strip_max > xdrop_best) xdrop_best = strip_max;
			else if (xdrop_best - strip_max >= xdrop) break;
		}
	}
	return best_i < 0 ? 0 : best_sc;
}

static Hsp antidiag_emit_traceback(const Sequence* q, Strand strand, int dna_len, const DpTarget& target,
	int best_sc, double evalue, int best_i, int best_j, int best_f, int j_beg, const AntidiagWorkspace& ws)
{
	Hsp out(true);
	out.swipe_target = target.target_idx; out.score = best_sc * config.cbs_matrix_scale;
	out.bit_score = score_matrix.bitscore(out.score); out.evalue = evalue;
	out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));
	int i = best_i, j = best_j, f = best_f, rem = best_sc;
	out.set_end(i + 1, j + 1, Frame(strand, f), dna_len);
	const int gap_open = score_matrix.gap_open(), gap_ext = score_matrix.gap_extend(), fs_pen = score_matrix.frame_shift();
	auto get_dir = [&](int r, int c, int fr) -> uint8_t {
		const int m_idx = (c - j_beg) >> 4, u = (c - j_beg) & 15;
		const AntidiagStripMeta& sm = ws.strips[(size_t)m_idx];
		const int s = r + u - sm.i_base;
		if (UNLIKELY(m_idx < 0 || (size_t)m_idx >= ws.strips.size() || s < 0 || s >= sm.s_end))
			throw runtime_error("Antidiag traceback bounds.");
		uint16_t d16;
		std::memcpy(&d16, ws.dir.data() + sm.dir_off + ((size_t)s * 16 + (size_t)u) * 2, sizeof(uint16_t));
		return (uint8_t)((d16 >> (fr * 5)) & 31);
	};
	while (rem > 0) {
		const uint8_t d = get_dir(i, j, f);
		if (!(d & 1)) {
			const Letter ql = letter_mask(q[f][i]), sl = letter_mask(target.seq[j]);
			const int m = score_matrix(sl, ql);
			out.push_match(ql, sl, m > 0);
			if (!(d & 2)) { rem -= m; --i; --j; }
			else if (!(d & 4)) { out.transcript.push_back(op_frameshift_forward); rem -= m - fs_pen; --i; --j; if (--f < 0) { f = 2; --i; } }
			else { out.transcript.push_back(op_frameshift_reverse); rem -= m - fs_pen; --i; --j; if (++f == 3) { f = 0; ++i; } }
		} else if (!(d & 4)) {
			int len = 1;
			while (get_dir(i, j - len, f) & 8) ++len;
			out.push_gap(op_deletion, len, target.seq.data() + j);
			rem += gap_open + len * gap_ext; j -= len;
		} else {
			int len = 1;
			while (get_dir(i - len, j, f) & 16) ++len;
			out.push_gap(op_insertion, len, nullptr);
			rem += gap_open + len * gap_ext; i -= len;
		}
	}
	if (UNLIKELY(rem != 0)) throw runtime_error("Antidiag traceback score mismatch.");
	out.set_begin(i + 1, j + 1, Frame(strand, f), dna_len);
	out.transcript.reverse();
	out.transcript.push_terminator();
	out.subject_source_range = out.subject_range;
	return out;
}

static list<Hsp> banded_3frame_antidiag(const TranslatedSequence& query, Strand strand,
	vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end, DpStat& stat, bool score_only)
{
	list<Hsp> out;
	if (target_begin == target_end) return out;
	Sequence q[3]; query.get_strand(strand, q);
	const int qlen = (int)q[0].length(), dna_len = (int)query.source().length();
	if (qlen <= 0 || dna_len < 3) return out;
	if (score_matrix.frame_shift() > 16384) {
		vector<DpTarget> ov32;
		return banded_3frame_swipe_targets<int32_t>(target_begin, target_end, score_only, query, strand, stat, false, ov32);
	}
	auto q_span = [&](const DpTarget& t) -> pair<int, int> {
		const int64_t jb = std::max<int64_t>(0, 1 - (int64_t)t.d_end);
		const int64_t je = std::min<int64_t>(t.seq.length(), (int64_t)qlen - t.d_begin);
		if (jb >= je || t.d_end <= t.d_begin) return { qlen, qlen };
		return { std::max(0, (int)jb + t.d_begin - 16), std::min(qlen, (int)je + t.d_end + 16) };
	};
	std::sort(target_begin, target_end, [&](const DpTarget& a, const DpTarget& b) { return q_span(a).first < q_span(b).first; });
	int min_slen = INT_MAX, max_slen = 0;
	for (auto it = target_begin; it < target_end; ++it) {
		const int sl = (int)it->seq.length();
		min_slen = std::min(min_slen, sl);
		max_slen = std::max(max_slen, sl);
	}
	int cutoff_sc = 1;
	if (min_slen <= max_slen) {
		int lo = 1, hi = 4096;
		while (lo < hi) {
			const int mid = (lo + hi) >> 1;
			const int sc_val = mid * config.cbs_matrix_scale;
			if (score_matrix.report_cutoff(sc_val, score_matrix.evalue(sc_val, qlen, (unsigned)min_slen)) ||
				score_matrix.report_cutoff(sc_val, score_matrix.evalue(sc_val, qlen, (unsigned)max_slen)))
				hi = mid;
			else
				lo = mid + 1;
		}
		cutoff_sc = lo;
	}
	const int xdrop = config.raw_swipe_xdrop;
	static thread_local AntidiagWorkspace ws;
	AntidiagCache& sc = ws.cache[strand == REVERSE ? 1 : 0];
	bool cache_hit = (sc.mat == score_matrix.matrix8());
	for (int f = 0; f < 3 && cache_hit; ++f)
		cache_hit = ((Loc)sc.q_copy[f].size() == q[f].length() && std::memcmp(sc.q_copy[f].data(), q[f].data(), (size_t)q[f].length() * sizeof(Letter)) == 0);
	if (!cache_hit) {
		if (sc.prof_buf.capacity() > (1u << 24)) vector<int8_t>().swap(sc.prof_buf);
		sc.mat = score_matrix.matrix8();
		for (int f = 0; f < 3; ++f) sc.q_copy[f].assign(q[f].data(), q[f].data() + q[f].length());
		sc.clusters.clear();
		sc.prof_used = 0;
		sc.hits.clear();
	} else if (score_only) {
		sc.hits.clear();
	}
	auto hit_less = [](const AntidiagPass1Hit& a, const AntidiagPass1Hit& b) {
		if (a.seq != b.seq) return (uintptr_t)a.seq < (uintptr_t)b.seq;
		if (a.d_begin != b.d_begin) return a.d_begin < b.d_begin;
		return a.d_end < b.d_end;
	};
	if (!score_only && sc.hits.size() > 1)
		std::sort(sc.hits.begin(), sc.hits.end(), hit_less);

	vector<DpTarget> ov32;
	for (auto c_beg = target_begin; c_beg < target_end;) {
		pair<int, int> sp = q_span(*c_beg);
		if (sp.first >= sp.second) { ++c_beg; continue; }
		const AntidiagCluster* hit = nullptr;
		for (size_t i = 0; i < sc.clusters.size(); ++i)
			if (sc.clusters[i].q0 <= sp.first && sp.second <= sc.clusters[i].q1) { hit = &sc.clusters[i]; break; }
		int q0 = sp.first, q1 = sp.second, pstride = 0;
		const int8_t* prof_ptr = nullptr;
		auto c_end = c_beg + 1;
		if (hit) {
			q0 = hit->q0; pstride = hit->pstride; prof_ptr = sc.prof_buf.data() + hit->prof_off;
			for (; c_end < target_end; ++c_end) {
				const pair<int, int> nsp = q_span(*c_end);
				if (nsp.first >= nsp.second || nsp.second > hit->q1) break;
			}
		} else {
			for (; c_end < target_end; ++c_end) {
				const pair<int, int> nsp = q_span(*c_end);
				if (nsp.first >= nsp.second || nsp.first > q1 + 256 || (std::max(q1, nsp.second) - q0 > 4096 && nsp.first > q0 + 64)) break;
				q1 = std::max(q1, nsp.second);
			}
			using Sv8 = ScoreVector<int8_t, 0>;
			pstride = ((q1 - q0 + 31) & ~31) + 96;
			const size_t prof_off = sc.prof_used, prof_sz = (size_t)AMINO_ACID_COUNT * 3 * (size_t)pstride;
			sc.prof_used += prof_sz;
			if (sc.prof_buf.size() < sc.prof_used) sc.prof_buf.resize(sc.prof_used);
			sc.clusters.push_back({ q0, q1, pstride, prof_off });
			int8_t* cl_prof = sc.prof_buf.data() + prof_off;
			std::memset(cl_prof, 0, prof_sz);
			prof_ptr = cl_prof;
			for (int f = 0; f < 3; ++f) {
				const int qfl = (int)q[f].length(), n = std::max(0, std::min(q1, qfl) - q0), full = n & ~31;
				const Letter* qp = q[f].data() + std::min(q0, qfl);
				alignas(32) Letter tail[32] = {}; std::copy(qp + full, qp + n, tail);
				const __m256i vtail = _mm256_load_si256((const __m256i*)tail);
				for (unsigned a = 0; a < AMINO_ACID_COUNT; ++a) {
					int8_t* dst = cl_prof + ((int)a * 3 + f) * pstride + 32;
					for (int i = 0; i < full; i += 32) Sv8(a, _mm256_loadu_si256((const __m256i*)(qp + i))).store(dst + i);
					if (full < n) Sv8(a, vtail).store(dst + full);
				}
			}
		}
		for (auto it = c_beg; it < c_end; ++it) {
			const int64_t jb = std::max<int64_t>(0, 1 - (int64_t)it->d_end);
			const int64_t je = std::min<int64_t>(it->seq.length(), (int64_t)qlen - it->d_begin);
			if (jb >= je || it->d_end <= it->d_begin) continue;
			const unsigned slen = (unsigned)it->seq.length();
			int t_cut = cutoff_sc;
			while (t_cut > 1 && score_matrix.report_cutoff((t_cut - 1) * config.cbs_matrix_scale, score_matrix.evalue((t_cut - 1) * config.cbs_matrix_scale, qlen, slen)))
				--t_cut;
			const int16_t init_best = (int16_t)(t_cut - 1);
			int best_i = 0, best_j = 0, best_f = 0;
			int best_sc = 0;
			if (score_only) {
				best_sc = antidiag_run_target_avx2<false>(q, prof_ptr, pstride, q0, *it, jb, je, qlen - 1, init_best, xdrop, ws, best_i, best_j, best_f);
			} else {
				const AntidiagPass1Hit key = { it->seq.data(), it->d_begin, it->d_end, 0, 0, 0 };
				const auto h_it = std::lower_bound(sc.hits.begin(), sc.hits.end(), key, hit_less);
				const AntidiagPass1Hit* h1 = (h_it != sc.hits.end() && h_it->seq == key.seq && h_it->d_begin == key.d_begin && h_it->d_end == key.d_end) ? &(*h_it) : nullptr;
				if (h1) {
					const int64_t je_crop = std::min<int64_t>(je, jb + (((h1->best_j - jb) | 15) + 1));
					const int i_crop = std::min(qlen - 1, (int)h1->best_i + 15);
					best_sc = antidiag_run_target_avx2<true>(q, prof_ptr, pstride, q0, *it, jb, je_crop, i_crop, (int16_t)(h1->best_sc - 1), 0, ws, best_i, best_j, best_f);
					if (UNLIKELY(best_sc != h1->best_sc))
						best_sc = antidiag_run_target_avx2<true>(q, prof_ptr, pstride, q0, *it, jb, je, qlen - 1, init_best, xdrop, ws, best_i, best_j, best_f);
				} else {
					best_sc = antidiag_run_target_avx2<true>(q, prof_ptr, pstride, q0, *it, jb, je, qlen - 1, init_best, xdrop, ws, best_i, best_j, best_f);
				}
			}
			if (UNLIKELY(best_sc >= SHRT_MAX)) {
				out.splice(out.end(), banded_3frame_swipe_targets<int32_t>(it, it + 1, score_only, query, strand, stat, false, ov32));
				continue;
			}
			if (best_sc <= 0) continue;
			const int score = best_sc * config.cbs_matrix_scale;
			const double evalue = score_matrix.evalue(score, qlen, slen);
			if (!score_matrix.report_cutoff(score, evalue)) continue;
			if (score_only) {
				sc.hits.push_back({ it->seq.data(), it->d_begin, it->d_end, best_i, best_j, (int16_t)best_sc });
				Hsp h(false);
				h.swipe_target = it->target_idx; h.score = score; h.bit_score = score_matrix.bitscore(score); h.evalue = evalue;
				h.query_range.end_ = std::min(best_i + 1, qlen);
				h.query_range.begin_ = std::max(h.query_range.end_ - (best_j + 1), 0);
				h.frame = strand == FORWARD ? 0 : 3;
				h.query_source_range = TranslatedPosition::absolute_interval(
					TranslatedPosition(h.query_range.begin_, Frame(h.frame)), TranslatedPosition(h.query_range.end_, Frame(h.frame)), dna_len);
				out.push_back(std::move(h));
			} else {
				out.push_back(antidiag_emit_traceback(q, strand, dna_len, *it, best_sc, evalue, best_i, best_j, best_f, (int)jb, ws));
			}
		}
		c_beg = c_end;
	}
	ws.trim(score_only);
	return out;
}
#endif

list<Hsp> banded_3frame_swipe(const TranslatedSequence &query, Strand strand, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end, DpStat &stat, bool score_only, bool parallel)
{
	bool pre_sorted = false;
#if defined(__AVX2__) && ARCH_AVX2_KERNELS
	if (config.antidiag_swipe && !parallel) {
		if (!score_only)
			return banded_3frame_antidiag(query, strand, target_begin, target_end, stat, false);
		const ptrdiff_t n_tgt = target_end - target_begin;
		const int qlen = (int)query.index(strand == FORWARD ? 0 : 3).length();
		if (n_tgt <= 8 || qlen <= 0)
			return banded_3frame_antidiag(query, strand, target_begin, target_end, stat, true);
		std::stable_sort(target_begin, target_end);
		pre_sorted = true;
		vector<DpTarget> anti_targets;
		auto swipe_end = target_begin;
		for (auto g = target_begin; g < target_end;) {
			const int m = (int)std::min<ptrdiff_t>(16, target_end - g);
			int band = 0, i1 = INT_MAX, q_lo = qlen, q_hi = 0;
			int64_t w_anti = 0;
			for (int k = 0; k < m; ++k) {
				const DpTarget& t = g[k];
				const int b = t.d_end - t.d_begin;
				band = std::max(band, b);
				i1 = std::min(i1, std::max(t.d_end - 1, 0));
				const int64_t jb = std::max<int64_t>(0, 1 - (int64_t)t.d_end);
				const int64_t je = std::min<int64_t>(t.seq.length(), (int64_t)qlen - t.d_begin);
				if (jb < je && b > 0) {
					w_anti += ((je - jb + 15) >> 4) * std::min(qlen, b + 15);
					q_lo = std::min(q_lo, std::max(0, (int)jb + t.d_begin - 16));
					q_hi = std::max(q_hi, std::min(qlen, (int)je + t.d_end + 16));
				}
			}
			int max_span = 0;
			for (int k = 0; k < m; ++k)
				max_span = std::max(max_span, (int)g[k].seq.length() - (i1 - (g[k].d_end - 1)));
			const int j_cols = std::max(0, std::min(max_span, qlen - (i1 + 1 - band)));
			const int64_t w_swipe = (int64_t)j_cols * std::min(qlen, band);
			if (q_hi > q_lo) w_anti += (int64_t)(q_hi - q_lo) * 2;
			if (m > 8 && w_swipe * ANTIDIAG_SWIPE_WEIGHT <= w_anti * 8) {
				if (swipe_end != g) std::move(g, g + m, swipe_end);
				swipe_end += m;
			} else {
				if (anti_targets.empty()) anti_targets.reserve((size_t)(target_end - g));
				anti_targets.insert(anti_targets.end(), std::make_move_iterator(g), std::make_move_iterator(g + m));
			}
			g += m;
		}
		if (swipe_end == target_begin)
			return banded_3frame_antidiag(query, strand, anti_targets.begin(), anti_targets.end(), stat, true);
		if (!anti_targets.empty()) {
			vector<DpTarget> ov16, ov32;
			list<Hsp> out = banded_3frame_swipe_targets<ScoreVector<int16_t, SHRT_MIN>>(target_begin, swipe_end, true, query, strand, stat, false, ov16);
			if (!ov16.empty())
				out.splice(out.end(), banded_3frame_swipe_targets<int32_t>(ov16.begin(), ov16.end(), true, query, strand, stat, false, ov32));
			out.splice(out.end(), banded_3frame_antidiag(query, strand, anti_targets.begin(), anti_targets.end(), stat, true));
			return out;
		}
	}
#endif
	vector<DpTarget> overflow16, overflow32;
#if defined (__SSE2__) | defined(__ARM_NEON)
	TaskTimer timer("Banded 3frame swipe (sort)", parallel ? 3 : UINT_MAX);
	if (!pre_sorted)
		std::stable_sort(target_begin, target_end);
	list<Hsp> out;
	if (parallel) {
		timer.go("Banded 3frame swipe (run)");
		vector<thread> threads;
		vector<list<Hsp>*> thread_out;
		vector<vector<DpTarget>> thread_overflow(config.threads_);
		atomic<size_t> next(0);
		for (int i = 0; i < config.threads_; ++i) {
			thread_out.push_back(new list<Hsp>);
			threads.emplace_back(banded_3frame_swipe_worker,
				target_begin,
				target_end,
				&next,
				score_only,
				&query,
				strand,
				thread_out.back(),
				&thread_overflow[i]);
		}
		for (auto &t : threads)
			t.join();
		timer.go("Banded 3frame swipe (merge)");
		for (list<Hsp>* l : thread_out) {
			out.splice(out.end(), *l);
			delete l;
		}
		overflow16.reserve(std::accumulate(thread_overflow.begin(), thread_overflow.end(), (size_t)0, [](size_t n, const vector<DpTarget> &v) { return n + v.size(); }));
		for (const vector<DpTarget> &v : thread_overflow)
			overflow16.insert(overflow16.end(), v.begin(), v.end());
	}
	else {
		if(score_only)
			out = banded_3frame_swipe_targets<ScoreVector<int16_t, SHRT_MIN>>(target_begin, target_end, score_only, query, strand, stat, false, overflow16);
		else
			out = banded_3frame_swipe_targets<int32_t>(target_begin, target_end, score_only, query, strand, stat, false, overflow16);
	}

	out.splice(out.end(), banded_3frame_swipe_targets<int32_t>(overflow16.begin(), overflow16.end(), score_only, query, strand, stat, false, overflow32));
	return out;
#else
	(void)pre_sorted;
	return banded_3frame_swipe_targets<int32_t>(target_begin, target_end, score_only, query, strand, stat, false, overflow32);
#endif
}

}

DISPATCH_7(std::list<Hsp>, banded_3frame_swipe, const TranslatedSequence&, query, Strand, strand, std::vector<DpTarget>::iterator, target_begin, std::vector<DpTarget>::iterator, target_end, DpStat&, stat, bool, score_only, bool, parallel)