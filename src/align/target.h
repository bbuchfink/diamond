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
#include <array>
#include <vector>
#include <list>
#include <memory>
#include "util/memory/memory_resource.h"
#include "util/geo/diagonal_segment.h"
#include "basic/const.h"
#include "util/hsp/approx_hsp.h"
#include "stats/hauser_correction.h"
#include "extend.h"
#include "util/data_structures/flat_array.h"
#include "stats/cbs.h"
#include "dp/flags.h"
#include "data/block/block.h"
#include "util/parallel/thread_pool.h"
#include "dp/score_profile.h"

struct SequenceSet;

namespace Search {

struct Hit;

}

namespace Extension {

struct SeedHit {
	int diag() const {
		return i - j;
	}
	bool operator<(const SeedHit& x) const {
		const int d1 = diag(), d2 = x.diag();
		return d1 < d2 || (d1 == d2 && j < x.j);
	}
	Interval query_range() const {
		return { i,i + 1 };
	}
	Interval target_range() const {
		return { j,j + 1 };
	}
	DiagonalSegment diag_segment() const {
		return DiagonalSegment(i, j, 1, score);
	}
	int i, j, score;
	unsigned frame;
};

struct WorkTarget {
	WorkTarget(BlockId block_id, const Sequence& seq, const Query& query, Loc max_target_len, Statistics& stats, std::pmr::memory_resource& pool);
	BlockId block_id;
	Sequence seq;
	std::array<int, MAX_CONTEXT> ungapped_score;
	std::array<std::list<ApproxHsp>, MAX_CONTEXT> hsp;
	std::unique_ptr<::Stats::TargetMatrix> matrix;
	LongScoreProfile<int16_t> profile, profile_rev;
	bool done;
};

WorkTarget ungapped_stage(std::vector<SeedHit>::iterator begin, std::vector<SeedHit>::iterator end, const Query& query, uint32_t block_id,
	Loc max_target_len, Statistics& stats, const Block& targets, const Mode mode, std::pmr::memory_resource& pool, const Search::Config& cfg);
std::vector<WorkTarget> ungapped_stage(const Query& query, FlatArray<SeedHit>::Iterator seed_hits, FlatArray<SeedHit>::Iterator seed_hits_end,
	std::vector<uint32_t>::const_iterator target_block_ids, DP::Flags flags, Statistics& stat, const Block& target_block, const Mode mode,
	std::pmr::memory_resource& pool, const Search::Config& cfg);

struct Target {

	Target(const BlockId block_id, const Sequence &seq, int ungapped_score, std::unique_ptr<::Stats::TargetMatrix>&& matrix):
		block_id(block_id),
		seq(seq),
		filter_score(0),
		filter_evalue(DBL_MAX),
		best_context(0),
		ungapped_score(ungapped_score),
		matrix(std::move(matrix)),
		done(false)
	{
	}

	void add_hit(Hsp&& hsp) {
		if(hsp.evalue < filter_evalue) {
			filter_evalue = hsp.evalue;
			filter_score = hsp.score;
			best_context = hsp.frame;
		}
		this->hsp[hsp.frame].push_back(std::move(hsp));
	}

	void add_hit(std::list<Hsp> &list, std::list<Hsp>::iterator it) {
		std::list<Hsp> &l = hsp[it->frame];
		l.splice(l.end(), list, it);
		if (l.back().score > filter_score) { // should be evalue
			filter_evalue = l.back().evalue;
			filter_score = l.back().score;
			best_context = l.back().frame;
		}
	}

	void add_hit(const ApproxHsp& h, Loc qlen) {
		hsp[h.frame].emplace_back(h, qlen, seq.length());
		if (h.evalue < filter_evalue) {
			filter_evalue = h.evalue;
			filter_score = h.score;
		}
		done = true;
	}

	static bool comp_evalue(const Target &t, const Target& u) {
		return t.filter_evalue < u.filter_evalue || (t.filter_evalue == u.filter_evalue && comp_score(t, u));
	}

	static bool comp_score(const Target& t, const Target& u) {
		return t.filter_score > u.filter_score || (t.filter_score == u.filter_score && t.block_id < u.block_id);
	}

	void inner_culling();
	void max_hsp_culling();

	BlockId block_id;
	Sequence seq;
	int filter_score;
	double filter_evalue;
	int best_context;
	int ungapped_score;
	std::array<std::list<Hsp>, MAX_CONTEXT> hsp;
	std::unique_ptr<::Stats::TargetMatrix> matrix;
	bool done;
};

using TargetList = std::pmr::vector<Target>;

struct TargetScore {
	uint32_t target;
	uint16_t score;
#ifdef EVAL_TARGET
	double evalue;
#endif
	bool operator<(const TargetScore& x) const {
#ifdef EVAL_TARGET
		return evalue < x.evalue || (evalue == x.evalue && target < x.target);
#else
		return score > x.score || (score == x.score && target < x.target);
#endif
	}
};

struct SeedHitList {
	FlatArray<SeedHit> seed_hits;
	std::vector<uint32_t> target_block_ids;
	std::vector<TargetScore> target_scores;
};

void culling(TargetList& targets, bool sort_only, const Search::Config& cfg);
void culling(std::vector<Match>& targets, const Search::Config& cfg);
bool append_hits(TargetList& targets, TargetList::iterator begin, TargetList::iterator end, const bool with_culling, const Search::Config& cfg);
std::pair<FlatArray<SeedHit>, std::vector<uint32_t>> gapped_filter(const Query& query, FlatArray<SeedHit>::Iterator seed_hits, FlatArray<SeedHit>::Iterator seed_hits_end,
	std::vector<uint32_t>::const_iterator target_block_ids, Statistics& stat, DP::Flags flags, const Search::Config &params);
TargetList align(std::vector<WorkTarget> &targets, const Query& query, DP::Flags flags, const HspValues hsp_values, const Mode mode, ThreadPool& tp, const Search::Config& cfg,
	Statistics &stat, std::pmr::memory_resource& pool);
std::vector<Match> align(TargetList &targets, const int64_t previous_matches, const Query& query, DP::Flags flags, const HspValues first_round, const bool first_round_culling,
	Statistics &stat, const Search::Config& cfg, std::pmr::memory_resource& pool);
TargetList full_db_align(const Query& query, DP::Flags flags, const HspValues hsp_values, Statistics &stat, const Block& target_block, std::pmr::memory_resource& pool);
void recompute_alt_hsps(std::vector<Match>::iterator begin, std::vector<Match>::iterator end, const Query& query, const HspValues v, Statistics& stats,
	std::pmr::memory_resource& pool);
void apply_filters(std::vector<Match>::iterator begin, std::vector<Match>::iterator end, const Query& query, const Search::Config& cfg);

std::vector<Match> extend(
	const Query& query,
	const Search::Config& cfg,
	Statistics &stat,
	DP::Flags flags,
	SeedHitList &l,
	std::pmr::memory_resource& pool);

}