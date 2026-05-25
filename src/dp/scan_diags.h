/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include "score_profile.h"

namespace DP {

void scan_diags128(const LongScoreProfile<int8_t>& qp, Sequence s, int d_begin, int j_begin, int j_end, int* out);
void scan_diags64(const LongScoreProfile<int8_t>& qp, Sequence s, int d_begin, int j_begin, int j_end, int* out);
void scan_diags(const LongScoreProfile<int8_t>& qp, Sequence s, int d_begin, int d_end, int j_begin, int j_end, int* out);
int diag_alignment(const int* s, int count);

}