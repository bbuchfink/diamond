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

#include <stdint.h>
#include "seed_array.h"
#include "enum_seeds.h"
#include "seed_array_impl.h"

using std::vector;
using std::fill;
using std::array;

namespace DISPATCH_ARCH {

template char* SeedArray<PackedLoc>::alloc_buffer(const SeedHistogram&, int);

template SeedArray<PackedLoc>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange&, int, char* buffer, const NoFilter*, const EnumCfg&);

}