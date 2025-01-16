/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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
#include "const.h"
#include "config.h"
#include "value.h"
#include "stats/score_matrix.h"

using PackedSeed = uint64_t;
#ifdef LONG_SEEDS
using SeedOffset = uint64_t;
using SeedPartition = uint32_t;
#else
using SeedOffset = uint32_t;
using SeedPartition = uint32_t;
#endif

static inline PackedSeed seedp_mask(int seedp_bits) {
	return ((PackedSeed)1 << seedp_bits) - 1;
}

static inline PackedSeed seedp_count(int seedp_bits) {
	return (PackedSeed)1 << seedp_bits;
}

static inline SeedPartition seed_partition(PackedSeed s, PackedSeed seedp_mask)
{
	return (SeedPartition)(s & seedp_mask);
}

static inline SeedOffset seed_partition_offset(PackedSeed s, PackedSeed seedp_bits)
{
	return (SeedOffset)(s >> seedp_bits);
}

struct Seed
{
	Letter& operator[](unsigned i)
	{
		return data_[i];
	}
	friend std::ostream& operator<<(std::ostream &str, const Seed &s)
	{
		for (unsigned i = 0; i < config.seed_weight; ++i)
			str << value_traits.alphabet[(size_t)s.data_[i]];
		return str;
	}
	int score(const Seed &rhs) const
	{
		int s = 0;
		for (unsigned i = 0; i < config.seed_weight; ++i)
			s += score_matrix(data_[i], rhs.data_[i]);
		return s;
	}
	uint64_t packed() const
	{
		uint64_t s = 0;
		for (unsigned i = 0; i < config.seed_weight; ++i) {
			s *= 20;
			s += data_[i];			
		}
		return s;
	}
	void enum_neighborhood(int treshold, std::vector<Seed> &out);
private:
	void enum_neighborhood(unsigned pos, int treshold, std::vector<Seed>& out, int score);
	Letter data_[Const::max_seed_weight];
};