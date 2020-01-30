/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

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

#include <list>
#include <atomic>
#include "../dp.h"
#include "../score_vector_int16.h"

using std::list;
using std::atomic;

namespace DP { namespace BandedSwipe { namespace DISPATCH_ARCH {

template<typename _sv, typename _traceback>
list<Hsp> swipe(
	const sequence &query,
	Frame frame,
	vector<DpTarget>::const_iterator subject_begin,
	vector<DpTarget>::const_iterator subject_end,
	const int8_t *composition_bias,
	int score_cutoff,
	vector<DpTarget> &overflow);

template<typename _sv>
list<Hsp> swipe_targets(const sequence &query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	int score_cutoff,
	vector<DpTarget> &overflow)
{
	list<Hsp> out;
	for (vector<DpTarget>::const_iterator i = begin; i < end; i += ScoreTraits<_sv>::CHANNELS) {
		if (flags & TRACEBACK)
			out.splice(out.end(), swipe<_sv, Traceback>(query, frame, i, i + std::min(vector<DpTarget>::const_iterator::difference_type(ScoreTraits<_sv>::CHANNELS), end - i), composition_bias, score_cutoff, overflow));
		else
			out.splice(out.end(), swipe<_sv, ScoreOnly>(query, frame, i, i + std::min(vector<DpTarget>::const_iterator::difference_type(ScoreTraits<_sv>::CHANNELS), end - i), composition_bias, score_cutoff, overflow));
	}
	return out;
}

void swipe_worker(const sequence *query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	atomic<size_t> *next,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	int score_cutoff,
	list<Hsp> *out,
	vector<DpTarget> *overflow)
{
	/*DpStat stat;
	size_t pos;
	vector<DpTarget> of;
	while (begin + (pos = next->fetch_add(config.swipe_chunk_size)) < end)
#ifdef __SSE2__
		out->splice(out->end(), banded_3frame_swipe_targets<score_vector<int16_t>>(begin + pos, min(begin + pos + config.swipe_chunk_size, end), score_only, *query, strand, stat, true, of));
#else
		out->splice(out->end(), banded_3frame_swipe_targets<int32_t>(begin + pos, min(begin + pos + config.swipe_chunk_size, end), score_only, *query, strand, stat, true, of));
#endif
	*overflow = std::move(of);*/
}


list<Hsp> swipe(const sequence &query, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end, Frame frame, const Bias_correction *composition_bias, int flags, int score_cutoff)
{
	vector<DpTarget> overflow16, overflow32;
#ifdef __SSE2__
	list<Hsp> out;
	std::stable_sort(target_begin, target_end);
	out = swipe_targets<score_vector<int16_t>>(query, target_begin, target_end, frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, score_cutoff, overflow16);
	if (!overflow16.empty())
		out.splice(out.end(), swipe_targets<int32_t>(query, overflow16.begin(), overflow16.end(), frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, score_cutoff, overflow32));
	return out;
#else
	return swipe_targets<int32_t>(query, target_begin, target_end, frame, composition_bias ? composition_bias->int8.data() : nullptr, flags, overflow32);
#endif
}
		
}}}