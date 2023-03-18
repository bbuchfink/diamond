/****
DIAMOND protein aligner
Copyright (C) 2019-2023 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#include "stage1.h"
#include "../util/simd/dispatch.h"

namespace Search {
	
DISPATCH_3V(run_stage1, JoinIterator<SeedLoc>&, it, Search::WorkSet*, work_set, const Search::Config*, cfg)

}