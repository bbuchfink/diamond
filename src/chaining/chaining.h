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
#include <utility>
#include <list>
#include "util/hsp/approx_hsp.h"
#include "dp/dp.h"

namespace Chaining {

std::pair<int, std::list<ApproxHsp>> run(Sequence query, Sequence subject, std::vector<DiagonalSegment>::const_iterator begin, std::vector<DiagonalSegment>::const_iterator end, bool log, unsigned frame);
std::list<Hsp> run(Sequence query, const std::vector<DpTarget>& targets);
ApproxHsp hamming_ext(std::vector<DiagonalSegment>::iterator begin, std::vector<DiagonalSegment>::iterator end, Loc qlen, Loc tlen, bool use_cov_filter);

}