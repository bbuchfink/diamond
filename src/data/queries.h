/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include <mutex>
#include <memory>
#include <vector>
#include "../util/io/output_file.h"
#include "block/block.h"

extern std::mutex query_aligned_mtx;
extern std::vector<bool> query_aligned;

void write_unaligned(const Block& query, OutputFile *file);
void write_aligned(const Block& query, OutputFile *file);

struct HashedSeedSet;
struct SeedSet;
extern std::unique_ptr<HashedSeedSet> query_seeds_hashed;
extern std::unique_ptr<SeedSet> query_seeds_bitset;