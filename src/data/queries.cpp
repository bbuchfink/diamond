/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include "queries.h"
#include "../util/sequence/sequence.h"
#include "../basic/config.h"
#include "seed_set.h"

using std::unique_ptr;

std::vector<bool> query_aligned;
std::mutex query_aligned_mtx;
unique_ptr<HashedSeedSet> query_seeds_hashed;
unique_ptr<SeedSet> query_seeds_bitset;

void write_unaligned(const Block& query, OutputFile *file)
{
	const size_t n = query.ids().size();
	for (size_t i = 0; i < n; ++i) {
		if (!query_aligned[i]) {
			Util::Seq::format(align_mode.query_translated ? query.source_seqs()[i] : query.seqs()[i],
				query.ids()[i],
				query.qual().empty() ? nullptr : query.qual()[i],
				*file,
				config.unfmt,
				input_value_traits);
		}
	}
}

void write_aligned(const Block& query, OutputFile *file)
{
	const size_t n = query.ids().size();
	for (size_t i = 0; i < n; ++i) {
		if (query_aligned[i]) {
			Util::Seq::format(align_mode.query_translated ? query.source_seqs()[i] : query.seqs()[i],
				query.ids()[i],
				query.qual().empty() ? nullptr : query.qual()[i],
				*file,
				config.alfmt,
				input_value_traits);
		}
	}
}