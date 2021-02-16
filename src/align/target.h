/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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
#include <array>
#include <set>
#include <vector>
#include <stdint.h>
#include <list>
#include <mutex>
#include <float.h>
#include "../search/trace_pt_buffer.h"
#include "../basic/diagonal_segment.h"
#include "../basic/const.h"
#include "../dp/hsp_traits.h"
#include "../stats/hauser_correction.h"
#include "extend.h"
#include "../util/data_structures/flat_array.h"
#include "../basic/parameters.h"
#include "../stats/cbs.h"

namespace Extension {

extern std::vector<int16_t*> target_matrices;
extern std::mutex target_matrices_lock;
extern std::atomic<size_t> target_matrix_count;

struct SeedHit {
	int diag() const {
		return i - j;
	}
	bool operator<(const SeedHit& x) const {
		const int d1 = diag(), d2 = x.diag();
		return d1 < d2 || (d1 == d2 && j < x.j);
	}
	int i, j, score;
	unsigned frame;
};

struct WorkTarget {
	WorkTarget(size_t block_id, const Sequence& seq, int query_len, const Stats::Composition& query_comp, const int16_t** query_matrix) :
		block_id(block_id),
		seq(seq)
	{
		ungapped_score.fill(0);
		if (config.comp_based_stats == Stats::CBS::HAUSER_AND_AVG_MATRIX_ADJUST) {
			const int l = (int)seq.length();
			const auto c = Stats::composition(seq);
			auto r = Stats::s_TestToApplyREAdjustmentConditional(query_len, l, query_comp.data(), c.data(), score_matrix.background_freqs());
			if (r == Stats::eCompoScaleOldMatrix)
				return;
			if (*query_matrix == nullptr) {
				*query_matrix = Stats::make_16bit_matrix(Stats::CompositionMatrixAdjust(query_len, query_len, query_comp.data(), query_comp.data(), Stats::CBS::AVG_MATRIX_SCALE, score_matrix.ideal_lambda(), score_matrix.joint_probs(), score_matrix.background_freqs()));
				++target_matrix_count;
			}
			if (target_matrices[block_id] == nullptr) {
				int16_t* target_matrix = Stats::make_16bit_matrix(Stats::CompositionMatrixAdjust(l, l, c.data(), c.data(), Stats::CBS::AVG_MATRIX_SCALE, score_matrix.ideal_lambda(), score_matrix.joint_probs(), score_matrix.background_freqs()));
				bool del = false;
				{
					std::lock_guard<std::mutex> lock(target_matrices_lock);
					if (target_matrices[block_id] == nullptr)
						target_matrices[block_id] = target_matrix;
					else del = true;
				}
				if (del)
					delete[] target_matrix;
				++target_matrix_count;
			}
			matrix = Stats::TargetMatrix(*query_matrix, target_matrices[block_id]);
		}
		else
			matrix = Stats::TargetMatrix(query_comp, query_len, seq);
	}
	bool adjusted_matrix() const {
		return !matrix.scores.empty();
	}
	size_t block_id;
	Sequence seq;
	std::array<int, MAX_CONTEXT> ungapped_score;
	std::array<std::list<Hsp_traits>, MAX_CONTEXT> hsp;
	Stats::TargetMatrix matrix;
};

std::vector<WorkTarget> ungapped_stage(const Sequence* query_seq, const Bias_correction* query_cb, const Stats::Composition& query_comp, FlatArray<SeedHit>& seed_hits, const std::vector<uint32_t>& target_block_ids, int flags, Statistics& stat);

struct Target {

	Target(size_t block_id, const Sequence &seq, int ungapped_score, const Stats::TargetMatrix& matrix):
		block_id(block_id),
		seq(seq),
		filter_score(0),
		filter_evalue(DBL_MAX),
		ungapped_score(ungapped_score),
		matrix(matrix)
	{}

	void add_hit(std::list<Hsp> &list, std::list<Hsp>::iterator it) {
		std::list<Hsp> &l = hsp[it->frame];
		l.splice(l.end(), list, it);
		filter_evalue = std::min(filter_evalue, l.back().evalue);
		filter_score = std::max(filter_score, l.back().score);
	}

	static bool comp_evalue(const Target &t, const Target& u) {
		return t.filter_evalue < u.filter_evalue || (t.filter_evalue == u.filter_evalue && comp_score(t, u));
	}

	static bool comp_score(const Target& t, const Target& u) {
		return t.filter_score > u.filter_score || (t.filter_score == u.filter_score && t.block_id < u.block_id);
	}

	bool adjusted_matrix() const {
		return !matrix.scores.empty();
	}

	void apply_filters(int source_query_len, const char *query_title, const Sequence& query_seq);
	void inner_culling(int source_query_len);
	void max_hsp_culling();

	size_t block_id;
	Sequence seq;
	int filter_score;
	double filter_evalue;
	int ungapped_score;
	std::array<std::list<Hsp>, MAX_CONTEXT> hsp;
	Stats::TargetMatrix matrix;
};

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

void load_hits(hit* begin, hit* end, FlatArray<SeedHit> &hits, std::vector<uint32_t> &target_block_ids, std::vector<TargetScore> &target_scores, unsigned query_len);
void culling(std::vector<Target>& targets, int source_query_len, const char* query_title, const Sequence& query_seq, size_t min_keep);
bool append_hits(std::vector<Target>& targets, std::vector<Target>::const_iterator begin, std::vector<Target>::const_iterator end, size_t chunk_size, int source_query_len, const char* query_title, const Sequence& query_seq);
std::vector<WorkTarget> gapped_filter(const Sequence *query, const Bias_correction* query_cbs, std::vector<WorkTarget>& targets, Statistics &stat);
void gapped_filter(const Sequence* query, const Bias_correction* query_cbs, FlatArray<SeedHit> &seed_hits, std::vector<uint32_t> &target_block_ids, Statistics& stat, int flags, const Parameters &params);
std::vector<Target> align(const std::vector<WorkTarget> &targets, const Sequence *query_seq, const Bias_correction *query_cb, int source_query_len, int flags, Statistics &stat);
std::vector<Match> align(std::vector<Target> &targets, const Sequence *query_seq, const Bias_correction *query_cb, int source_query_len, int flags, Statistics &stat);
std::vector<Target> full_db_align(const Sequence *query_seq, const Bias_correction *query_cb, int flags, Statistics &stat);

std::vector<Match> extend(
	size_t query_id,
	const Parameters &params,
	const Metadata &metadata,
	Statistics &stat,
	int flags,
	FlatArray<SeedHit>& seed_hits,
	std::vector<uint32_t>& target_block_ids,
	const std::vector<TargetScore>& target_scores);

struct Memory {
	Memory(size_t query_count);
	int& low_score(size_t query_id);
	int& mid_score(size_t query_id);
	int& min_score(size_t query_id, size_t i);
	size_t count(size_t query_id) const;
	void update(size_t query_id, std::vector<Target>::const_iterator begin, std::vector<Target>::const_iterator end);
	void update_failed_count(size_t query_id, size_t failed_count, int ranking_low_score);
	int ranking_low_score(size_t query_id) const {
		return ranking_low_score_[query_id];
	}
	size_t ranking_failed_count(size_t query_id) const {
		return ranking_failed_count_[query_id];
	}
private:
	const size_t N;
	std::vector<int> scores_, count_, ranking_low_score_, ranking_failed_count_;
};

extern Memory* memory;

}
