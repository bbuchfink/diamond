/****
DIAMOND protein aligner
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
#include <vector>
#include "hamming/hit_field.h"
#include "basic/statistics.h"
#include "util/algo/pattern_matcher.h"
#include "util/memory/alignment.h"
#include "run/config.h"
#include "hit_buffer.h"
#include "hit.h"
#include "data/flags.h"
#include "kmer_ranking.h"
#include "util/algo/join_result.h"
#include "basic/reduction.h"
#include "util/data_structures/deque.h"
#include "data/block/block.h"

// #define UNGAPPED_SPOUGE

struct SensitivityTraits {
	const bool     support_query_indexed;
	const bool     motif_masking;
	const double   freq_sd;
	const unsigned min_identities;
	const double   ungapped_evalue;
	const double   ungapped_evalue_short;
	const double   gapped_filter_evalue;
	const unsigned index_chunks;
	const unsigned query_bins;
	const char*    contiguous_seed;
	const double   seed_cut;
	const double   default_block_size;
	const Reduction& reduction;
	const int      minimizer_window;
#if WITH_DNA
  const double   chain_fraction_align;
  const int       min_chain_score;
  const double max_overlap_extension;
#endif
  const int      sketch_size;
};

struct HashedSeedSet;

namespace Search {

struct Context {
	const PatternMatcher previous_matcher, current_matcher;
	const int short_query_ungapped_cutoff;
	KmerRanking* kmer_ranking;
	const PackedSeed seedp_mask;
};

extern Reduction steinegger12;
extern Reduction murphy10;
extern Reduction no_reduction;
extern Reduction dna;
extern const std::map<Sensitivity, SensitivityTraits> sensitivity_traits;
extern const std::map<Sensitivity, std::vector<std::string>> shape_codes;
extern const std::map<Sensitivity, std::vector<Round>> iterated_sens;

void search_shape(unsigned sid, int query_block, unsigned query_iteration, char* query_buffer, char* ref_buffer, Config& cfg, const HashedSeedSet* target_seeds);
bool use_single_indexed(double coverage, size_t query_letters, size_t ref_letters);
void setup_search(Sensitivity sens, Search::Config& cfg);
MaskingAlgo soft_masking_algo(const SensitivityTraits& traits);
int seedp_bits(int shape_weight, int threads, int index_chunks);

}

namespace Search {

using Container = vector<std::array<char, 48>, Util::Memory::AlignmentAllocator<std::array<char, 48>, 16>>;

struct WorkSet {
	WorkSet(const Context& context, const Search::Config& cfg, unsigned shape_id, HitBuffer::Writer* out, AsyncWriter<Hit, Search::Config::RankingBuffer::EXPONENT>* global_ranking_buffer, KmerRanking *kmer_ranking):
		context(context),
		cfg(cfg),
		shape_id(shape_id),
		out(out),
		global_ranking_buffer(global_ranking_buffer),
		kmer_ranking(kmer_ranking)
	{}
	Context context;
	const Search::Config& cfg;
	unsigned shape_id;
	Statistics stats;
	HitBuffer::Writer* out;
	AsyncWriter<Hit, Search::Config::RankingBuffer::EXPONENT>* global_ranking_buffer;
#ifndef __APPLE__
	Container vq, vs;
#endif
	HitField hits;
	KmerRanking* kmer_ranking;
};

void run_stage1(JoinIterator<PackedLoc>& it, WorkSet* work_set, const Search::Config* cfg);
void run_stage1(JoinIterator<PackedLocId>& it, WorkSet* work_set, const Search::Config* cfg);
bool keep_target_id(const Search::Config& cfg);

}

extern const double SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE;
