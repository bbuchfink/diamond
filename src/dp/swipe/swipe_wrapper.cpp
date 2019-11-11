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
#include "../dp.h"
#include "../score_vector_int16.h"

using std::list;

namespace DP { namespace BandedSwipe { namespace DISPATCH_ARCH {

template<typename _sv, typename _traceback>
list<Hsp> swipe(
	const sequence &query,
	Frame frame,
	vector<DpTarget>::const_iterator subject_begin,
	vector<DpTarget>::const_iterator subject_end,
	const typename ScoreTraits<_sv>::Score *composition_bias,
	vector<DpTarget> &overflow);

template<typename _sv>
list<Hsp> swipe_targets(const sequence &query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	Frame frame,
	const typename ScoreTraits<_sv>::Score *composition_bias,
	int flags,
	vector<DpTarget> &overflow)
{
	list<Hsp> out;
	for (vector<DpTarget>::const_iterator i = begin; i < end; i += ScoreTraits<_sv>::CHANNELS) {
		if (flags & TRACEBACK)
			out.splice(out.end(), swipe<_sv, Traceback>(query, frame, i, i + std::min(vector<DpTarget>::const_iterator::difference_type(ScoreTraits<_sv>::CHANNELS), end - i), composition_bias, overflow));
		else
			out.splice(out.end(), swipe<_sv, ScoreOnly>(query, frame, i, i + std::min(vector<DpTarget>::const_iterator::difference_type(ScoreTraits<_sv>::CHANNELS), end - i), composition_bias, overflow));
	}
	return out;
}

list<Hsp> swipe(const sequence &query, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end, Frame frame, const Bias_correction &composition_bias, int flags)
{
	vector<DpTarget> overflow16, overflow32;
#ifdef __SSE2__
	list<Hsp> out;
	std::stable_sort(target_begin, target_end);
	out = swipe_targets<score_vector<int16_t>>(query, target_begin, target_end, frame, composition_bias.int16.data(), flags, overflow16);
	return out;
#endif
}
		
}}}