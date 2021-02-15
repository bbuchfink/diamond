/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include <vector>
#include "../score_vector.h"
#include "../score_vector_int8.h"
#include "../score_vector_int16.h"
#include "swipe.h"
#include "../../basic/sequence.h"
#include "target_iterator.h"
#include "../../util/data_structures/mem_buffer.h"

using std::vector;
using std::list;
using std::max;
using namespace DISPATCH_ARCH;

namespace DP { namespace Swipe {
namespace DISPATCH_ARCH {

template<typename _sv>
struct Matrix
{
	struct ColumnIterator
	{
		ColumnIterator(_sv* hgap_front, _sv* score_front) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_;
		}
		inline _sv hgap() const
		{
			return *hgap_ptr_;
		}
		inline _sv diag() const
		{
			return *score_ptr_;
		}
		inline void set_hgap(const _sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const _sv& x)
		{
			*score_ptr_ = x;
		}
		std::nullptr_t trace_mask() {
			return nullptr;
		}
		_sv *hgap_ptr_, *score_ptr_;
	};
	Matrix(int rows, int)
	{
		hgap_.resize(rows);
		score_.resize(rows + 1);
		std::fill(hgap_.begin(), hgap_.end(), ScoreTraits<_sv>::zero());
		std::fill(score_.begin(), score_.end(), ScoreTraits<_sv>::zero());
	}
	inline ColumnIterator begin(int)
	{
		return ColumnIterator(hgap_.begin(), score_.begin());
	}
	void set_zero(int c)
	{
		const int l = (int)hgap_.size();
		for (int i = 0; i < l; ++i) {
			set_channel(hgap_[i], c, ScoreTraits<_sv>::zero_score());
			set_channel(score_[i], c, ScoreTraits<_sv>::zero_score());
		}
		set_channel(score_[l], c, ScoreTraits<_sv>::zero_score());
	}
	constexpr int cols() const {
		return 1;
	}
private:
#ifdef __APPLE__
	MemBuffer<_sv> hgap_, score_;
#else
	static thread_local MemBuffer<_sv> hgap_, score_;
#endif
};

#ifndef __APPLE__
template<typename _sv> thread_local MemBuffer<_sv> Matrix<_sv>::hgap_;
template<typename _sv> thread_local MemBuffer<_sv> Matrix<_sv>::score_;
#endif

template<typename _sv>
struct TracebackVectorMatrix
{
	typedef typename ScoreTraits<_sv>::TraceMask TraceMask;
	typedef void* Stat;
	struct ColumnIterator
	{
		ColumnIterator(_sv* hgap_front, _sv* score_front, TraceMask* trace_mask_front) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front),
			trace_mask_ptr_(trace_mask_front)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_; ++trace_mask_ptr_;
		}
		inline _sv hgap() const
		{
			return *hgap_ptr_;
		}
		inline _sv diag() const
		{
			return *score_ptr_;
		}
		inline TraceMask* trace_mask() {
			return trace_mask_ptr_;
		}
		inline void set_hgap(const _sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const _sv& x)
		{
			*score_ptr_ = x;
		}
		std::nullptr_t stat() {
			return nullptr;
		}
		std::nullptr_t hstat() {
			return nullptr;
		}
		void set_hstat(std::nullptr_t) {}
		inline void set_zero() {}
		_sv *hgap_ptr_, *score_ptr_;
		TraceMask* trace_mask_ptr_;
	};

	struct TracebackIterator
	{
		TracebackIterator(const TraceMask *mask, const TraceMask* mask_begin, const TraceMask* mask_end, int rows, int i, int j, size_t channel) :
			rows_(rows),
			mask_(mask),
			mask_begin_(mask_begin),
			mask_end_(mask_end),
			channel_mask_vgap(TraceMask::vmask(channel)),
			channel_mask_hgap(TraceMask::hmask(channel)),
			i(i),
			j(j)
		{
			assert(i >= 0 && j >= 0);
		}
		void wrap_mask() {
			if (mask_ < mask_begin_)
				mask_ = mask_end_ - (mask_begin_ - mask_);
		}
		TraceMask mask() const {
			return *mask_;
		}
		void walk_diagonal()
		{
			mask_ -= rows_ + 1;
			wrap_mask();
			--i;
			--j;
			assert(i >= -1 && j >= -1);
		}
		pair<Edit_operation, int> walk_gap()
		{
			if (mask_->gap & channel_mask_vgap) {
				int l = 0;
				do {
					++l;
					--i;
					--mask_;
				} while (((mask_->open & channel_mask_vgap) == 0) && (i > 0));
				return std::make_pair(op_insertion, l);
			}
			else {
				int l = 0;
				do {
					++l;
					--j;
					mask_ -= rows_;
					wrap_mask();
				} while (((mask_->open & channel_mask_hgap) == 0) && (j > 0));
				return std::make_pair(op_deletion, l);
			}
		}
		const int rows_;
		const TraceMask* mask_, *mask_begin_, *mask_end_;
		const decltype(TraceMask::gap) channel_mask_vgap, channel_mask_hgap;
		int i, j;
	};

	TracebackIterator traceback(int col, int i, int j, size_t channel) const
	{
		return TracebackIterator(&trace_mask_[col*rows_ + i], trace_mask_.begin(), trace_mask_.end(), rows_, i, j, channel);
	}

	TracebackVectorMatrix(int rows, int cols) :
		rows_(rows),
		cols_(cols)
	{
		hgap_.resize(rows);
		score_.resize(rows + 1);
		trace_mask_.resize(cols * rows);
		std::fill(hgap_.begin(), hgap_.end(), _sv());
		std::fill(score_.begin(), score_.end(), _sv());

	}

	inline ColumnIterator begin(int col)
	{
		return ColumnIterator(hgap_.begin(), score_.begin(), &trace_mask_[col*rows_]);
	}

	void set_zero(int c)
	{
		const int l = (int)hgap_.size();
		for (int i = 0; i < l; ++i) {
			set_channel(hgap_[i], c, ScoreTraits<_sv>::zero_score());
			set_channel(score_[i], c, ScoreTraits<_sv>::zero_score());
		}
		set_channel(score_[l], c, ScoreTraits<_sv>::zero_score());
	}

	int cols() const {
		return cols_;
	}

#ifdef __APPLE__
	MemBuffer<_sv> hgap_, score_;
#else
	static thread_local MemBuffer<_sv> hgap_, score_;
#endif
	MemBuffer<TraceMask> trace_mask_;
private:
	int rows_, cols_;
};

#ifndef __APPLE__
template<typename _sv> thread_local MemBuffer<_sv> TracebackVectorMatrix<_sv>::hgap_;
template<typename _sv> thread_local MemBuffer<_sv> TracebackVectorMatrix<_sv>::score_;
#endif

template<typename _sv, typename _traceback>
struct MatrixTraits
{};

template<typename _sv>
struct MatrixTraits<_sv, VectorTraceback>
{
	typedef TracebackVectorMatrix<_sv> Type;
	typedef RowCounter<_sv> MyRowCounter;
};

template<typename _sv>
struct MatrixTraits<_sv, ScoreOnly>
{
	typedef Matrix<_sv> Type;
	typedef DummyRowCounter MyRowCounter;
};

template<typename _sv, typename _cbs>
Hsp traceback(const sequence& query, Frame frame, _cbs bias_correction, const Matrix<_sv>& dp, const DpTarget& target, typename ScoreTraits<_sv>::Score max_score, double evalue, int max_col, int max_i, int max_j, int channel)
{
	Hsp out;
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<_sv>::int_score(max_score) * config.cbs_matrix_scale;
	out.evalue = evalue;
	out.frame = frame.index();
	return out;
}

template<typename _sv, typename _cbs>
Hsp traceback(const sequence &query, Frame frame, _cbs bias_correction, const TracebackVectorMatrix<_sv> &dp, const DpTarget &target, typename ScoreTraits<_sv>::Score max_score, double evalue, int max_col, int max_i, int max_j, int channel)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	typedef typename ScoreTraits<_sv>::TraceMask TraceMask;
	const auto channel_mask = TraceMask::vmask(channel) | TraceMask::hmask(channel);
	typename TracebackVectorMatrix<_sv>::TracebackIterator it(dp.traceback(max_col, max_i, max_j, channel));
	Hsp out;
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<_sv>::int_score(max_score) * config.cbs_matrix_scale;
	out.evalue = evalue;
	out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));

	out.frame = frame.index();
	out.query_range.end_ = it.i + 1;
	out.subject_range.end_ = it.j + 1;
	const int end_score = out.score;
	int score = 0;

	while (it.i >= 0 && it.j >= 0 && score < end_score) {
		if ((it.mask().gap & channel_mask) == 0) {
			const Letter q = query[it.i], s = target.seq[it.j];
			const int m = score_matrix(q, s);
			const int m2 = add_cbs_scalar(m, bias_correction[it.i]);
			score += m2;
			out.push_match(q, s, m > (Score)0);
			it.walk_diagonal();
		}
		else {
			const pair<Edit_operation, int> g(it.walk_gap());
			out.push_gap(g.first, g.second, target.seq.data() + it.j + g.second);
			score -= score_matrix.gap_open() + g.second * score_matrix.gap_extend();
		}
	}

	if (score != end_score)
		throw std::runtime_error("Traceback error.");

	out.query_range.begin_ = it.i + 1;
	out.subject_range.begin_ = it.j + 1;
	out.transcript.reverse();
	out.transcript.push_terminator();
	return out;
}

template<typename _sv, typename _traceback, typename _cbs>
list<Hsp> swipe(const sequence& query, Frame frame, DynamicIterator<DpTarget>& target_it, _cbs composition_bias, vector<DpTarget>& overflow, Statistics &stats)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	typedef typename MatrixTraits<_sv, _traceback>::Type Matrix;
	constexpr int CHANNELS = ScoreTraits<_sv>::CHANNELS;

	int max_col[CHANNELS], max_i[CHANNELS], max_j[CHANNELS];
	const int qlen = (int)query.length();

	if (qlen > MatrixTraits<_sv, _traceback>::MyRowCounter::MAX_LEN)
		throw std::runtime_error("Query length exceeds row counter maximum.");

	const _sv open_penalty(static_cast<Score>(score_matrix.gap_open() + score_matrix.gap_extend())),
		extend_penalty(static_cast<Score>(score_matrix.gap_extend()));
	//_sv best = _sv();
	Score best[CHANNELS];
	std::fill(best, best + CHANNELS, ScoreTraits<_sv>::zero_score());
	SwipeProfile<_sv> profile;
	AsyncTargetBuffer<Score> targets(target_it);
	Matrix dp(qlen, targets.max_len());
	CBSBuffer<_sv, _cbs> cbs_buf(composition_bias, qlen, 0);
	list<Hsp> out;
	int col = 0;
	
	while (targets.active.size() > 0) {
		typename Matrix::ColumnIterator it(dp.begin(col));
		typename MatrixTraits<_sv, _traceback>::MyRowCounter row_counter(0);
		_sv vgap, hgap, last, col_best;
		vgap = hgap = last = col_best = _sv();
		profile.set(targets.seq_vector());
#ifdef DP_STAT
		stats.inc(Statistics::GROSS_DP_CELLS, uint64_t(qlen) * CHANNELS);
#endif
		for (int i = 0; i < qlen; ++i) {
			hgap = it.hgap();
			const _sv next = swipe_cell_update<_sv>(it.diag(), profile.get(query[i]), cbs_buf(i), extend_penalty, open_penalty, hgap, vgap, col_best, nullptr, nullptr, nullptr, it.trace_mask(), row_counter);
			it.set_hgap(hgap);
			it.set_score(last);
			last = next;
			++it;
		}
		it.set_score(last);
		//best = max(best, col_best);
		
		Score col_best_[CHANNELS], i_max[CHANNELS];
		store_sv(col_best, col_best_);
		row_counter.store(i_max);
		for (int i = 0; i < targets.active.size();) {
			int c = targets.active[i];
			if (col_best_[c] > best[c]) {
				best[c] = col_best_[c];
				max_col[c] = col;
				max_i[c] = ScoreTraits<_sv>::int_score(i_max[c]);
				max_j[c] = targets.pos[c];
			}
			bool reinit = false;
			if (col_best_[c] == ScoreTraits<_sv>::max_score()) {
				overflow.push_back(targets.dp_targets[c]);
				reinit = true;
			} else if (!targets.inc(c)) {
				const int s = ScoreTraits<_sv>::int_score(best[c]) * config.cbs_matrix_scale;
				const double evalue = score_matrix.evalue(s, qlen, (unsigned)targets.dp_targets[c].seq.length());
				if (score_matrix.report_cutoff(s, evalue))
					out.push_back(traceback<_sv>(query, frame, composition_bias, dp, targets.dp_targets[c], best[c], evalue, max_col[c], max_i[c], max_j[c], c));
				reinit = true;				
			}
			if (reinit) {
				if (targets.init_target(i, c)) {
					dp.set_zero(c);
					//set_channel(best, c, ScoreTraits<_sv>::zero_score());
					best[c] = ScoreTraits<_sv>::zero_score();
				}
				else
					continue;
			}
			++i;
		}
		col = (col + 1) % dp.cols();
	}

	return out;
}

#ifdef __SSE4_1__
template list<Hsp> swipe<score_vector<int8_t>, VectorTraceback, const int8_t*>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, const int8_t*, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<score_vector<int8_t>, ScoreOnly, const int8_t*>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, const int8_t*, vector<DpTarget>&, Statistics&);
#endif
#ifdef __SSE2__
template list<Hsp> swipe<score_vector<int16_t>, VectorTraceback, const int8_t*>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, const int8_t*, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<score_vector<int16_t>, ScoreOnly, const int8_t*>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, const int8_t*, vector<DpTarget>&, Statistics&);
#endif
template list<Hsp> swipe<int32_t, VectorTraceback, const int8_t*>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, const int8_t*, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<int32_t, ScoreOnly, const int8_t*>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, const int8_t*, vector<DpTarget>&, Statistics&);

#ifdef __SSE4_1__
template list<Hsp> swipe<score_vector<int8_t>, VectorTraceback, NoCBS>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, NoCBS, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<score_vector<int8_t>, ScoreOnly, NoCBS>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, NoCBS, vector<DpTarget>&, Statistics&);
#endif
#ifdef __SSE2__
template list<Hsp> swipe<score_vector<int16_t>, VectorTraceback, NoCBS>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, NoCBS, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<score_vector<int16_t>, ScoreOnly, NoCBS>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, NoCBS, vector<DpTarget>&, Statistics&);
#endif
template list<Hsp> swipe<int32_t, VectorTraceback, NoCBS>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, NoCBS, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<int32_t, ScoreOnly, NoCBS>(const sequence&, Frame, DynamicIterator<DpTarget>& target_it, NoCBS, vector<DpTarget>&, Statistics&);

}}}