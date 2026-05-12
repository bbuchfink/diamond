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
#include "basic/value.h"

namespace Geo {

inline Loc i(Loc j, Loc d) {
	return d + j;
}

inline Loc j(Loc i, Loc d) {
	return i - d;
}

inline Loc diag_sub_matrix(Loc d, Loc i0, Loc j0) {
	return d + j0 - i0;
}

inline Loc rev_diag(Loc d, Loc qlen, Loc tlen) {
	return -d + qlen - tlen;
}

inline Loc min_diag(Loc qlen, Loc tlen) {
	return -(tlen - 1);
}

inline Loc max_diag(Loc qlen, Loc tlen) {
	return qlen - 1;
}

inline Loc clip_diag(Loc d, Loc qlen, Loc tlen) {
	return std::min(std::max(d, min_diag(qlen, tlen)), max_diag(qlen, tlen));
}

inline void assert_diag_bounds(Loc d, Loc qlen, Loc tlen) {
	assert(d >= min_diag(qlen, tlen) && d <= max_diag(qlen, tlen));
}

}