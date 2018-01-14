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

#include "dp.h"

int sw_3frame(const TranslatedSequence &query, Strand strand, const sequence &subject, int gap_open, int gap_extend, int frame_shift, Hsp &out)
{
	using std::max;

	Matrix<int> score, hgap, vgap;
	score.init(int(query.source().length()-2), (int)subject.length());
	hgap.init(int(query.source().length()-2), (int)subject.length());
	vgap.init(int(query.source().length()-2), (int)subject.length());

	int max_score = 0, max_i = -1, max_j = -1;

	for (int i = 0; i < query.source().length()-2; ++i)
		for (int j = 0; j < subject.length(); ++j) {
			int s = 0;
			const int m = score_matrix(query(i, strand), subject[j]);
			s = max(s, ((i >= 3 && j >= 1) ? score[i - 3][j - 1] : 0) + m);
			s = max(s, (j >= 1 ? hgap[i][j - 1] : 0));
			s = max(s, (i >= 3 ? vgap[i - 3][j] : 0));
			if (i >= 4 && j >= 1)
				s = max(s, score[i - 4][j - 1] + m - frame_shift);
			if (i >= 2 && j >= 1)
				s = max(s, score[i - 2][j - 1] + m - frame_shift);
			score[i][j] = s;
			if (s > max_score) {
				max_score = s;
				max_i = i;
				max_j = j;
			}
			vgap[i][j] = max(i >= 3 ? vgap[i - 3][j] - gap_extend : 0, s - gap_open - gap_extend);
			hgap[i][j] = max(j >= 1 ? hgap[i][j - 1] - gap_extend : 0, s - gap_open - gap_extend);
		}

	int i = max_i;
	int j = max_j;
	int s;

	cout << '\t';
	for (int j = 0; j < subject.length(); ++j) {
		cout << j << '\t';
	}
	cout << endl;
	for (int i = 0; i < query.source().length() - 2; ++i) {
		cout << i << '\t';
		for (int j = 0; j < subject.length(); ++j) {
			cout << score[i][j] << '\t';
		}
		cout << endl;
	}
	return max_score;

	out.query_range.end_ = TranslatedPosition(i, strand).translated + 1;
	out.subject_range.end_ = j + 1;

	while (i >= 0 && j >= 0 && (s = score[i][j]) > 0) {
		const int m = score_matrix(query(i, strand), subject[j]);
		if (i >= 3 && j >= 1 && s == score[i - 3][j - 1] + m) {
			if (query(i, strand) == subject[j])
				out.transcript.push_back(op_match, 1u);
			else
				out.transcript.push_back(op_substitution, subject[j]);
			i -= 3;
			--j;
		}
		else if (j >= 1 && s == hgap[i][j - 1]) {
			out.transcript.push_back(op_deletion, subject[j]);
			--j;
		}
		else if (i >= 3 && s == vgap[i - 3][j]) {
			out.transcript.push_back(op_insertion, 1u);
			i -= 3;
		}
		else if (i >= 4 && j >= 1 && s == score[i - 4][j - 1] + m - frame_shift) {
			if (query(i, strand) == subject[j])
				out.transcript.push_back(op_match, 1u);
			else
				out.transcript.push_back(op_substitution, subject[j]);
			out.transcript.push_back(op_frameshift_forward);
			i -= 4;
			--j;
		}
		else if (i >= 2 && j >= 1 && s == score[i - 2][j - 1] + m - frame_shift) {
			if (query(i, strand) == subject[j])
				out.transcript.push_back(op_match, 1u);
			else
				out.transcript.push_back(op_substitution, subject[j]);
			out.transcript.push_back(op_frameshift_reverse);
			i -= 2;
			--j;
		}
		else
			throw std::runtime_error("Traceback error");
		++out.length;
	}

	out.query_range.begin_ = TranslatedPosition(i + 3, strand);
	out.frame = TranslatedPosition(i + 3, strand).frame.index();
	out.subject_range.begin_ = j + 1;
	out.transcript.reverse();
	out.transcript.push_terminator();

	return max_score;
}