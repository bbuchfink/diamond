/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include "../basic/shape.h"
#include "../util/data_structures/double_array.h"
#include "../data/flags.h"
#include "../run/config.h"
#include "../data/seed_histogram.h"

namespace Search {

struct SeedStats {
	SeedStats():
		good_seed_positions(0),
		low_complexity_seeds(0)
	{}
	size_t good_seed_positions, low_complexity_seeds;
};

bool seed_is_complex(const Letter* seq, const Shape& shape, const double cut);
bool seed_is_complex_unreduced(Letter* seq, const Shape& shape, const double cut, const bool mask_seeds, SeedStats& stats);
void mask_seeds(const Shape& shape, const SeedPartitionRange& range, DoubleArray<SeedLoc>* query_seed_hits, DoubleArray<SeedLoc>* ref_seed_hits, Search::Config& cfg);

}