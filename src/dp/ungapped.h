/****
DIAMOND protein aligner
Copyright (C) 2013-2022 Max Planck Society for the Advancement of Science e.V.
						Benjamin Buchfink
						Eberhard Karls Universitaet Tuebingen

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

#pragma once
#include "../basic/value.h"
#include "../util/geo/diagonal_segment.h"
#include "../stats/hauser_correction.h"

//int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned seed_len, unsigned &delta, unsigned &len);
//int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned &delta, unsigned &len);
int xdrop_ungapped_right(const Letter *query, const Letter *subject, int &len);
int ungapped_window(const Letter* query, const Letter* subject, int window);
DiagonalSegment xdrop_ungapped(const Sequence &query, const Bias_correction &query_bc, const Sequence &subject, int qa, int sa);
DiagonalSegment xdrop_ungapped(const Sequence& query, const int8_t* query_cbs, const Sequence& subject, int qa, int sa, bool count_identities);
DiagonalSegment xdrop_ungapped(const Sequence& query, const Sequence& subject, const DiagonalSegment& anchor);
int score_range(Sequence query, Sequence subject, int i, int j, int j_end);
template<typename Cbs>
DiagonalSegment score_range_s(Sequence query, Cbs query_cbs, Sequence subject, int i, int j, int j_end);
Score self_score(const Sequence& seq);
Hsp trivial(Sequence query, Sequence target, const int8_t* query_cbs);
Anchor make_clipped_anchor(const Anchor& anchor, Sequence query, const int8_t* query_cbs, Sequence target);
Anchor make_null_anchor(const Anchor& anchor);