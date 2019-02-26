/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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
#include "swipe_matrix.h"
#include "swipe.h"
#include "target_iterator.h"

#ifdef __SSE2__

template<typename _score>
void banded_swipe(const sequence &query, vector<DpTarget>::iterator subject_begin, vector<DpTarget>::iterator subject_end)
{
	/*typedef score_vector<_score> sv;

	assert(subject_end - subject_begin <= score_traits<_score>::channels);
	const int qlen = (int)query.length();

	int i0 = INT_MAX, i1 = INT_MIN, cols = 0;
	for (vector<DpTarget>::const_iterator j = subject_begin; j < subject_end; ++j) {
		int i2 = std::max(j->d_end - 1, 0);
		int j0 = i2 - (j->d_end - 1),
			j1 = std::min(qlen - 1 - j->d_begin, (int)(j->seq.length() - 1)) + 1;
		cols = std::max(cols, j1 - j0);
		i1 = std::max(i1, i2);
		i0 = std::min(i0, i2 + 1 - (j->d_end - j->d_begin));
	}
	const int band = i1 + 1 - i0;
	assert(band > 0);

	//BandedSwipeMatrix<_score> dp(band);
	BandedSwipeTracebackMatrix<_score> dp(band, cols);

	const sv open_penalty(static_cast<char>(score_matrix.gap_open() + score_matrix.gap_extend())),
		extend_penalty(static_cast<char>(score_matrix.gap_extend()));
	sv best;
	best.zero();
	SwipeProfile<_score> profile;
	TargetIterator<score_traits<_score>::channels> targets(subject_begin, subject_end, Banded());

	int j = 0;
	while (targets.active.size() > 0) {
		const int i0_ = std::max(i0, 0), i1_ = std::min(i1, qlen - 1);
		if (i0_ > i1_)
			break;
		typename BandedSwipeTracebackMatrix<_score>::ColumnIterator it(dp.begin(i0_ - i0, j));
		if (i0_ - i0 > 0)
			it.set_zero();
		sv vgap, hgap;
		vgap.zero();
		
		profile.set(targets.get());
		for (int i = i0_; i <= i1_; ++i) {
			hgap = it.hgap();
			const sv next = cell_update<_score>(it.diag(), profile.get(query[i]), extend_penalty, open_penalty, hgap, vgap, best);
			it.set_hgap(hgap);
			it.set_score(next);
			++it;
		}

		cells += targets.active.size()*(i1_ - i0_ + 1);
		for (int i = 0; i < targets.active.size(); ++i) {
			int channel = targets.active[i];
			if (!targets.inc(channel))
				targets.active.erase(i);
		}
		++i0;
		++i1;
		++j;
	}
	for (int i = 0; i < targets.n_targets; ++i)
		subject_begin[i].score = best[i];*/
}

#endif

void banded_swipe(const sequence &query, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end)
{
#ifdef __SSE2__
	banded_swipe<int16_t>(query, target_begin, target_end);
#endif
}