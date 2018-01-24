/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include <vector>
#include "../score_vector.h"
#include "swipe.h"
#include "swipe_matrix.h"
#include "../../basic/sequence.h"
#include "target_iterator.h"

// #define SW_ENABLE_DEBUG

using std::vector;
using std::pair;

#ifdef __SSE2__

template<typename _sv> TLS_PTR vector<_sv>* SwipeMatrix<_sv>::hgap_ptr;
template<typename _sv> TLS_PTR vector<_sv>* SwipeMatrix<_sv>::score_ptr;

template<typename _sv>
void swipe(const sequence &query, vector<DpTarget>::const_iterator subject_begin, vector<DpTarget>::const_iterator subject_end, vector<int>::iterator out)
{
#ifdef SW_ENABLE_DEBUG
	static int v[1024][1024];
#endif

	const int qlen = (int)query.length();
	SwipeMatrix<_sv> dp(qlen);

	const _sv open_penalty(static_cast<char>(score_matrix.gap_open() + score_matrix.gap_extend())),
		extend_penalty(static_cast<char>(score_matrix.gap_extend())),
		vbias(score_matrix.bias());
	_sv best;
	SwipeProfile<_sv> profile;
	TargetIterator<_sv::CHANNELS> targets(subject_begin, subject_end);

	while (targets.active.size() > 0) {
		typename SwipeMatrix<_sv>::ColumnIterator it(dp.begin());
		_sv vgap, hgap, last;
		//profile.set(targets.get2<_score>());
		for (int i = 0; i < qlen; ++i) {
			hgap = it.hgap();
			const _sv next = cell_update<_sv>(it.diag(), profile.get(query[i]), extend_penalty, open_penalty, hgap, vgap, best, vbias);
			it.set_hgap(hgap);
			it.set_score(last);
			last = next;
#ifdef SW_ENABLE_DEBUG
			v[targets.pos[0]][i] = next[0];
#endif
			++it;
		}
		it.set_score(last);
		
		for (int i = 0; i < targets.active.size();) {
			int j = targets.active[i];
			if (!targets.inc(j)) {
				out[targets.target[j]] = best[j];
				if (targets.init_target(i, j)) {
					dp.set_zero(j);
					best.set(j, 0);
				}
				else
					continue;
			}
			++i;
		}
	}

#ifdef SW_ENABLE_DEBUG
	for (unsigned j = 0; j < qlen; ++j) {
		for (unsigned i = 0; i < subject_begin[0].length(); ++i)
			printf("%4i", v[i][j]);
		printf("\n");
	}
	printf("\n");
#endif
}

#endif

void swipe(const sequence &query, vector<DpTarget>::const_iterator subject_begin, vector<DpTarget>::const_iterator subject_end, vector<int>::iterator out)
{
#ifdef __SSE2__
	swipe<score_vector<uint8_t> >(query, subject_begin, subject_end, out);
#endif
}