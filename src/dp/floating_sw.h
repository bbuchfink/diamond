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

#ifndef FLOATING_SW_H_
#define FLOATING_SW_H_

#include "../basic/match.h"
#include "../align/align.h"

struct Score_only { };
struct Traceback { };

template<typename _score, typename _traceback, typename _score_correction>
void floating_sw(const Letter *query,
	const Letter *subject,
	Hsp_data &segment,
	int band,
	_score xdrop,
	_score gap_open,
	_score gap_extend,
	uint64_t &cell_updates,
	unsigned query_anchor,
	unsigned subject_anchor,
	int min_j,
	const _score_correction &score_correction,
	const _traceback& = Score_only(),
	const _score& = int());

#endif /* FLOATING_SW_H_ */
