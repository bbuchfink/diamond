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

#include "../dp.h"

namespace DP { namespace BandedSwipe { namespace DISPATCH_ARCH {

template<typename _sv, typename _traceback> void swipe(const sequence &query, Frame frame, vector<DpTarget>::const_iterator subject_begin, vector<DpTarget>::const_iterator subject_end);

template<typename _sv>
void swipe_targets(const sequence &query,
	vector<DpTarget>::iterator begin,
	vector<DpTarget>::iterator end,
	Frame frame,
	int flags)
{
	for (vector<DpTarget>::iterator i = begin; i < end; i += ScoreTraits<_sv>::CHANNELS) {
		if (flags & TRACEBACK)
			swipe<_sv, Traceback>(query, frame, i, i + std::min(vector<DpTarget>::iterator::difference_type(ScoreTraits<_sv>::CHANNELS), end - i));
		else
			swipe<_sv, ScoreOnly>(query, frame, i, i + std::min(vector<DpTarget>::iterator::difference_type(ScoreTraits<_sv>::CHANNELS), end - i));
	}
}

void swipe(const sequence &query, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end, Frame frame, int flags)
{
#ifdef __SSE2__
	std::stable_sort(target_begin, target_end);
	swipe_targets<score_vector<int16_t>>(query, target_begin, target_end, frame, flags);
#endif
}
		
}}}