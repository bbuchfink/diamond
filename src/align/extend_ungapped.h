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

#ifndef EXTEND_UNGAPPED_H_
#define EXTEND_UNGAPPED_H_

#include "../dp/dp.h"
#include "../data/reference.h"

inline Diagonal_segment ungapped_extension(unsigned subject, unsigned subject_pos, unsigned query_pos, const sequence &query)
{
	const Letter* s = ref_seqs::data_->data(ref_seqs::data_->position(subject, subject_pos)),
		*q = &query[query_pos];
	unsigned delta, len;
	int score = xdrop_ungapped(q, s, delta, len);
	return Diagonal_segment(query_pos - delta, subject_pos - delta, len, score);
}

inline Diagonal_segment ungapped_extension(unsigned subject_pos, unsigned query_pos, const sequence &query, const sequence &subject)
{
	const Letter* s = &subject[subject_pos],
		*q = &query[query_pos];
	unsigned delta, len;
	int score = xdrop_ungapped(q, s, delta, len);
	return Diagonal_segment(query_pos - delta, subject_pos - delta, len, score);
}

#endif