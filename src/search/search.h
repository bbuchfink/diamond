/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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
bool keep_target_id(const Search::Config& cfg);
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

}

extern const double SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE;