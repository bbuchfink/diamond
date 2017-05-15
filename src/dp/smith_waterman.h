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

#ifndef SSE_SW_H_
#define SSE_SW_H_

#include <algorithm>
#include <stddef.h>
#include "score_profile.h"
#include "dp_matrix.h"

#ifdef __SSE2__

template<typename _score>
inline score_vector<_score> cell_update(const score_vector<_score> &diagonal_cell,
						 const score_vector<_score> &scores,
						 const score_vector<_score> &gap_extension,
						 const score_vector<_score> &gap_open,
						 score_vector<_score> &horizontal_gap,
						 score_vector<_score> &vertical_gap,
						 score_vector<_score> &best,
						 const score_vector<_score> &vbias)
{
	score_vector<_score> current_cell = diagonal_cell + scores;
	current_cell.unbias(vbias);
	current_cell.max(vertical_gap).max(horizontal_gap);
	best.max(current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -=  gap_extension;
	score_vector<_score> open = current_cell - gap_open;
	vertical_gap.max(open);
	horizontal_gap.max(open);
	return current_cell;
}

template<typename _score, typename _callback>
void smith_waterman(const sequence &query,
			const vector<sequence> &subjects,
			unsigned band,
			unsigned padding,
			int op,
			int ep,
			int filter_score,
			_callback &f,
			const _score&,
			Statistics &stats)
{
	#ifdef SW_ENABLE_DEBUG
	int v[1024][1024];
	#endif

	typedef score_vector<_score> sv;

	unsigned qlen ((unsigned)query.length());
	unsigned slen ((unsigned)subjects[0].length());
	DP_matrix<_score> dp (slen, qlen, band, padding);

	sv open_penalty (static_cast<char>(op));
	sv extend_penalty (static_cast<char>(ep));
	sv vbias (score_matrix.bias());
	sequence_stream dseq;
	score_profile<_score> profile;

	typename vector<sequence>::const_iterator subject_it (subjects.begin());
	while(subject_it < subjects.end()) {

		const unsigned n_subject (std::min((unsigned)score_traits<_score>::channels, (unsigned)(subjects.end() - subject_it)));
		typename vector<sequence>::const_iterator subject_end (subject_it + n_subject);
		sv best;
		dseq.reset();
		dp.clear();

		for(unsigned j=0;j<slen;++j) {
			typename DP_matrix<_score>::Column_iterator it (dp.begin(j));
			sv vgap, hgap, column_best;
			profile.set(dseq.get(subject_it, subject_end, j, _score()));

			while(!it.at_end()) {
				hgap = it.hgap();
				sv next = cell_update<_score>(it.diag(), profile.get(query[it.row_pos_]), extend_penalty, open_penalty, hgap, vgap, column_best, vbias);
				it.set_hgap(hgap);
				it.set_score(next);
				#ifdef SW_ENABLE_DEBUG
				v[j][it.row_pos_] = next[0];
				#endif
				++it;
			}
			best.max(column_best);
		}

		for(unsigned i=0;i<n_subject;++i)
			if(best[i] >= filter_score)
				f(i + unsigned(subject_it - subjects.begin()), *(subject_it + i), best[i]);
		subject_it += std::min((ptrdiff_t)score_traits<_score>::channels, subjects.end()-subject_it);
	}

	#ifdef SW_ENABLE_DEBUG
	for(unsigned j=0;j<qlen;++j) {
		for(unsigned i=0;i<subjects[0].length();++i)
			printf("%4i", v[i][j]);
		printf("\n");
	}
	printf("\n");
	#endif
}

#endif

#endif /* SSE_SW_H_ */
