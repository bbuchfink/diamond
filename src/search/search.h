/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
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
#include <stddef.h>
#include <vector>
#include "../util/data_structures/flat_array.h"
#include "../util/simd.h"
#include "../dp/ungapped.h"
#include "../basic/shape_config.h"
#include "../basic/statistics.h"
#include "../util/algo/pattern_matcher.h"
#include "../util/scores/cutoff_table.h"
#include "finger_print.h"
#include "../util/memory/alignment.h"
#include "../run/config.h"
#include "hit.h"
#include "../util/data_structures/writer.h"
#include "../data/seed_array.h"

// #define UNGAPPED_SPOUGE

struct SensitivityTraits {
	bool     support_query_indexed;
	double   freq_sd;
	unsigned min_identities;
	double   ungapped_evalue;
	double   ungapped_evalue_short;
	double   gapped_filter_evalue;
	unsigned index_chunks;
	unsigned query_bins;
};

extern const std::map<Sensitivity, SensitivityTraits> sensitivity_traits;
extern const std::map<Sensitivity, std::vector<std::string>> shape_codes;
extern const std::map<Sensitivity, std::vector<Sensitivity>> iterated_sens;

namespace Search {

struct Context {
	const PatternMatcher previous_matcher, current_matcher;
#ifdef UNGAPPED_SPOUGE
	const Util::Scores::CutoffTable2D cutoff_table;
#else
	const Util::Scores::CutoffTable cutoff_table, cutoff_table_short;
#endif
	const int short_query_ungapped_cutoff;
};

}

struct Stage1_hit
{
	Stage1_hit(unsigned q_ref, unsigned q_offset, unsigned s_ref, unsigned s_offset) :
		q(q_ref + q_offset),
		s(s_ref + s_offset)
	{}
	bool operator<(const Stage1_hit &rhs) const
	{
		return q < rhs.q;
	}
	struct Query
	{
		unsigned operator()(const Stage1_hit &x) const
		{
			return x.q;
		}
	};
	unsigned q, s;
};

struct HashedSeedSet;

void search_shape(unsigned sid, unsigned query_block, unsigned query_iteration, char *query_buffer, char *ref_buffer, Search::Config& cfg, const HashedSeedSet* target_seeds);
bool use_single_indexed(double coverage, size_t query_letters, size_t ref_letters);
void setup_search(Sensitivity sens, Search::Config& cfg);

namespace Search {

struct WorkSet {
	Context context;
	const Search::Config& cfg;
	unsigned shape_id;
	Statistics stats;
	Writer<Hit>* out;
#ifndef __APPLE__
	std::vector<FingerPrint, Util::Memory::AlignmentAllocator<FingerPrint, 16>> vq, vs;
#endif
	FlatArray<uint32_t> hits;
};

DECL_DISPATCH(void, stage1, (const SeedArray::Entry::Value* q, size_t nq, const SeedArray::Entry::Value* s, size_t ns, WorkSet& work_set))

}

extern const double SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE;