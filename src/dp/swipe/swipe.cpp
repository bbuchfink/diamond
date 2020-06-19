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
		_sv *hgap_ptr_, *score_ptr_;
	};
	Matrix(int rows)
	{
		hgap_.resize(rows);
		score_.resize(rows + 1);
		std::fill(hgap_.begin(), hgap_.end(), ScoreTraits<_sv>::zero());
		std::fill(score_.begin(), score_.end(), ScoreTraits<_sv>::zero());
	}
	inline ColumnIterator begin()
	{
		return ColumnIterator(&hgap_[0], &score_[0]);
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
private:
	static thread_local MemBuffer<_sv> hgap_, score_;
};

template<typename _sv> thread_local MemBuffer<_sv> Matrix<_sv>::hgap_;
template<typename _sv> thread_local MemBuffer<_sv> Matrix<_sv>::score_;

template<typename _sv, typename _cbs>
Hsp traceback(const sequence& query, Frame frame, _cbs bias_correction, const Matrix<_sv>& dp, const DpTarget& target, typename ScoreTraits<_sv>::Score max_score, int max_col, int channel)
{
	Hsp out;
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<_sv>::int_score(max_score);
	out.frame = frame.index();
	return out;
}

template<typename _sv, typename _traceback, typename _cbs>
//list<Hsp> swipe(const sequence &query, const sequence *subject_begin, const sequence *subject_end, int score_cutoff, vector<int> &overflow)
list<Hsp> swipe(const sequence& query, Frame frame, vector<DpTarget>::const_iterator subject_begin, vector<DpTarget>::const_iterator subject_end, _cbs composition_bias, int score_cutoff, vector<DpTarget>& overflow, Statistics &stats)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	constexpr int CHANNELS = ScoreTraits<_sv>::CHANNELS;

	int max_col[CHANNELS];
	const int qlen = (int)query.length();
	Matrix<_sv> dp(qlen);

	const _sv open_penalty(static_cast<Score>(score_matrix.gap_open() + score_matrix.gap_extend())),
		extend_penalty(static_cast<Score>(score_matrix.gap_extend()));
	_sv best = _sv();
	SwipeProfile<_sv> profile;
	TargetBuffer<Score> targets(subject_begin, subject_end);
	list<Hsp> out;

	while (targets.active.size() > 0) {
		typename Matrix<_sv>::ColumnIterator it(dp.begin());
		_sv vgap, hgap, last, col_best;
		vgap = hgap = last = col_best = _sv();
		profile.set(targets.seq_vector());
#ifdef DP_STAT
		stats.inc(Statistics::GROSS_DP_CELLS, uint64_t(qlen) * CHANNELS);
#endif
		for (int i = 0; i < qlen; ++i) {
			hgap = it.hgap();
			const _sv next = swipe_cell_update<_sv>(it.diag(), profile.get(query[i]), nullptr, extend_penalty, open_penalty, hgap, vgap, col_best);
			it.set_hgap(hgap);
			it.set_score(last);
			last = next;
			++it;
		}
		it.set_score(last);
		best = max(best, col_best);
		
		Score col_best_[CHANNELS];
		store_sv(col_best, col_best_);
		for (int i = 0; i < targets.active.size();) {
			int c = targets.active[i];
			if (col_best_[c] == ScoreTraits<_sv>::max_score()) {
				overflow.push_back(targets.dp_target(c));
				if (targets.init_target(i, c)) {
					dp.set_zero(c);
					set_channel(best, c, ScoreTraits<_sv>::zero_score());
				}
				else
					continue;
			} else if (!targets.inc(c)) {
				const int s = ScoreTraits<_sv>::int_score(extract_channel(best, c));
				if (s >= score_cutoff)
					out.push_back(traceback<_sv>(query, frame, composition_bias, dp, targets.dp_target(c), extract_channel(best, c), max_col[c], c));
				if (targets.init_target(i, c)) {
					dp.set_zero(c);
					set_channel(best, c, ScoreTraits<_sv>::zero_score());
				}
				else
					continue;
			}
			++i;
		}
	}

	return out;
}

/*list<Hsp> swipe(const sequence &query, const sequence *subject_begin, const sequence *subject_end, int score_cutoff)
{
	vector<int> overflow8, overflow16, overflow32;
#ifdef __SSE4_1__
	list<Hsp> out = swipe<score_vector<int8_t>>(query, subject_begin, subject_end, score_cutoff, overflow8);

	if (overflow8.empty())
		return out;

	vector<sequence> overflow_seq;
	overflow_seq.reserve(overflow8.size());
	for (int i : overflow8)
		overflow_seq.push_back(subject_begin[i]);
	list<Hsp> out16 = swipe<score_vector<int16_t>>(query, overflow_seq.data(), overflow_seq.data() + overflow_seq.size(), score_cutoff, overflow16);
	for (Hsp &hsp : out16)
		hsp.swipe_target = overflow8[hsp.swipe_target];
	out.splice(out.end(), out16);

	if (overflow16.empty())
		return out;

	overflow_seq.clear();
	for (int i : overflow16)
		overflow_seq.push_back(subject_begin[overflow8[i]]);
	list<Hsp> out32 = swipe<int32_t>(query, overflow_seq.data(), overflow_seq.data() + overflow_seq.size(), score_cutoff, overflow32);
	for (Hsp &hsp : out32)
		hsp.swipe_target = overflow8[overflow16[hsp.swipe_target]];
	out.splice(out.end(), out32);

	return out;
#else
	return {};
#endif
}*/

#ifdef __SSE4_1__
template list<Hsp> swipe<score_vector<int8_t>, Traceback, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<score_vector<int8_t>, ScoreOnly, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);
#endif
#ifdef __SSE2__
template list<Hsp> swipe<score_vector<int16_t>, Traceback, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<score_vector<int16_t>, ScoreOnly, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);
#endif
template list<Hsp> swipe<int32_t, Traceback, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<int32_t, ScoreOnly, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);

}}}