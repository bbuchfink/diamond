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

#include <algorithm>
#include <utility>
#include <math.h>
#include <mutex>
#include "extend.h"
#include "../data/queries.h"
#include "../basic/config.h"
#include "../dp/comp_based_stats.h"
#include "target.h"
#include "../dp/dp.h"
#include "../util/log_stream.h"
#include "../data/reference.h"
#include "../util/system.h"
#include "culling.h"
#include "../util/util.h"

using std::vector;
using std::list;
using std::array;
using std::pair;

namespace Extension {

constexpr uint64_t MAX_CHUNK_SIZE = 400, MIN_CHUNK_SIZE = 128;

size_t ranking_chunk_size(size_t target_count) {
	if (config.no_ranking)
		return target_count;
	if (config.ext_chunk_size > 0)
		return config.ext_chunk_size;
	const double letters = ref_seqs::get().letters(), default_letters = config.sensitivity >= Sensitivity::VERY_SENSITIVE ? 400 * 1e6 : 2 * 1e9;
	const uint64_t block_mult = std::max(uint64_t(letters / default_letters), (uint64_t)1);
	if (config.toppercent < 100.0)
		return MIN_CHUNK_SIZE;
	return std::max(MIN_CHUNK_SIZE, std::min(make_multiple(config.max_alignments, (uint64_t)32), MAX_CHUNK_SIZE));
}

size_t chunk_size_multiplier(const FlatArray<SeedHit>& seed_hits, int query_len) {
	return seed_hits.size() * query_len / seed_hits.data_size() < config.seedhit_density ? config.chunk_size_multiplier : 1;
}

struct TargetScore {
	uint32_t target;
	uint16_t score;
	bool operator<(const TargetScore& x) const {
		return score > x.score || (score == x.score && target < x.target);
	}
};

void load_hits(hit* begin, hit* end, FlatArray<SeedHit> &hits, vector<uint32_t> &target_block_ids, vector<TargetScore> &target_scores) {
	hits.clear();
	hits.reserve(end - begin);
	target_block_ids.clear();
	target_scores.clear();
	if (begin >= end)
		return;
	std::sort(begin, end, hit::CmpSubject());
	const size_t total_subjects = ref_seqs::get().get_length();
	uint32_t target = UINT32_MAX;
	uint16_t score = 0;
	if (std::log2(total_subjects) * (end - begin) < total_subjects / 10 ) {
		for (hit* i = begin; i < end; ++i) {
			std::pair<size_t, size_t> l = ref_seqs::data_->local_position((uint64_t)i->subject_);
			const uint32_t t = (uint32_t)l.first;
			if (t != target) {
				if (target != UINT32_MAX) {
					target_scores.push_back({ uint32_t(target_block_ids.size() - 1), score });
					score = 0;
				}
				hits.next();
				target = t;
				target_block_ids.push_back(target);
			}
			hits.push_back({ (int)i->seed_offset_, (int)l.second, i->query_ % align_mode.query_contexts });
			score = std::max(score, i->score_);
		}
	} else {
		typename vector<size_t>::const_iterator limit_begin = ref_seqs::get().limits_begin(), it = limit_begin;
		for (const hit* i = begin; i < end; ++i) {
			const size_t subject_offset = (uint64_t)i->subject_;
			while (*it <= subject_offset) ++it;
			uint32_t t = (uint32_t)(it - limit_begin) - 1;
			if (t != target) {
				if (target != UINT32_MAX) {
					target_scores.push_back({ uint32_t(target_block_ids.size() - 1), score });
					score = 0;
				}
				hits.next();
				target_block_ids.push_back(t);
				target = t;
			}
			hits.push_back({ (int)i->seed_offset_, (int)(subject_offset - *(it - 1)), i->query_ % align_mode.query_contexts });
			score = std::max(score, i->score_);
		}
	}
	if (target != UINT32_MAX)
		target_scores.push_back({ uint32_t(target_block_ids.size() - 1), score });
}

vector<Target> extend(const Parameters& params,
	size_t query_id,
	const sequence *query_seq,
	int source_query_len,
	const Bias_correction *query_cb,
	FlatArray<SeedHit> &seed_hits,
	vector<uint32_t> &target_block_ids,
	const Metadata& metadata,
	Statistics& stat,
	int flags)
{
	stat.inc(Statistics::TARGET_HITS1, target_block_ids.size());
	task_timer timer(flags & DP::PARALLEL ? config.target_parallel_verbosity : UINT_MAX);
	if (config.gapped_filter_evalue > 0.0) {
		timer.go("Computing gapped filter");
		gapped_filter(query_seq, query_cb, seed_hits, target_block_ids, stat, flags, params);
		if ((flags & DP::PARALLEL) == 0)
			stat.inc(Statistics::TIME_GAPPED_FILTER, timer.microseconds());
	}
	stat.inc(Statistics::TARGET_HITS2, target_block_ids.size());

	timer.go("Computing chaining");
	vector<WorkTarget> targets = ungapped_stage(query_seq, query_cb, seed_hits, target_block_ids.data(), flags);
	if ((flags & DP::PARALLEL) == 0)
		stat.inc(Statistics::TIME_CHAINING, timer.microseconds());

	return align(targets, query_seq, query_cb, source_query_len, flags, stat);
}

vector<Match> extend(const Parameters &params, size_t query_id, hit* begin, hit* end, const Metadata &metadata, Statistics &stat, int flags) {
	const unsigned contexts = align_mode.query_contexts;
	vector<sequence> query_seq;
	vector<Bias_correction> query_cb;
	const char* query_title = query_ids::get()[query_id];

	if (config.log_query || flags & DP::PARALLEL)
		log_stream << "Query=" << query_title << " Hits=" << end - begin << endl;

	for (unsigned i = 0; i < contexts; ++i)
		query_seq.push_back(query_seqs::get()[query_id*contexts + i]);

	task_timer timer(flags & DP::PARALLEL ? config.target_parallel_verbosity : UINT_MAX);
	if (config.comp_based_stats == 1) {
		timer.go("Computing CBS");
		for (unsigned i = 0; i < contexts; ++i)
			query_cb.emplace_back(query_seq[i]);
		timer.finish();
	}

	const int source_query_len = align_mode.query_translated ? (int)query_source_seqs::get()[query_id].length() : (int)query_seqs::get()[query_id].length();

	timer.go("Loading seed hits");
	TLS_FIX_S390X FlatArray<SeedHit> seed_hits, seed_hits_chunk;
	thread_local vector<uint32_t> target_block_ids, target_block_ids_chunk;
	thread_local vector<TargetScore> target_scores;
	load_hits(begin, end, seed_hits, target_block_ids, target_scores);
	stat.inc(Statistics::TARGET_HITS0, target_block_ids.size());
	stat.inc(Statistics::TIME_LOAD_HIT_TARGETS, timer.microseconds());
	timer.finish();
	const size_t target_count = target_block_ids.size();
	
	const size_t chunk_size = ranking_chunk_size(target_count);

	if (chunk_size < target_count) {
		timer.go("Sorting targets by score");
		std::sort(target_scores.begin(), target_scores.end());
		stat.inc(Statistics::TIME_SORT_TARGETS_BY_SCORE, timer.microseconds());
		timer.finish();
	}
		
	const int relaxed_cutoff = score_matrix.rawscore(config.min_bit_score == 0.0
		? score_matrix.bitscore(config.max_evalue * config.relaxed_evalue_factor, (unsigned)query_seq[0].length())
		: config.min_bit_score);
	vector<TargetScore>::const_iterator i0 = target_scores.cbegin(), i1 = std::min(i0 + chunk_size, target_scores.cend());
	if (config.toppercent == 100.0)
		while (i1 < target_scores.cend() && i1->score >= relaxed_cutoff && (i1 - i0) < config.max_alignments) ++i1;
	const int low_score = config.query_memory ? memory->low_score(query_id) : 0;
	const size_t previous_count = config.query_memory ? memory->count(query_id) : 0;
	bool first_round_traceback = config.min_id > 0 || config.query_cover > 0 || config.subject_cover > 0;
	//size_t multiplier = 1;
	int tail_score = 0;
	if (first_round_traceback)
		flags |= DP::TRACEBACK;

	vector<Target> aligned_targets;
	while(i0 < target_scores.cend()) {
		seed_hits_chunk.clear();
		target_block_ids_chunk.clear();
		const size_t current_chunk_size = (size_t)(i1 - i0);
		const bool multi_chunk = current_chunk_size < target_scores.size();
		if (config.query_memory && memory->ranking_failed_count(query_id) >= chunk_size && memory->ranking_low_score(query_id) >= i0->score)
			break;

		if (multi_chunk) {
			for (vector<TargetScore>::const_iterator j = i0; j < i1; ++j) {
				target_block_ids_chunk.push_back(target_block_ids[j->target]);
				seed_hits_chunk.push_back(seed_hits.begin(j->target), seed_hits.end(j->target));
			}
		}
		else {
			target_block_ids_chunk = TLS_FIX_S390X_MOVE(target_block_ids);
			seed_hits_chunk = TLS_FIX_S390X_MOVE(seed_hits);
		}

		//multiplier = std::max(multiplier, chunk_size_multiplier(seed_hits_chunk, (int)query_seq.front().length()));

		vector<Target> v = extend(params, query_id, query_seq.data(), source_query_len, query_cb.data(), seed_hits_chunk, target_block_ids_chunk, metadata, stat, flags);
		const size_t n = v.size();
		stat.inc(Statistics::TARGET_HITS3, v.size());
		bool new_hits = false;
		if (multi_chunk)
			new_hits = append_hits(aligned_targets, v.begin(), v.end(), chunk_size, source_query_len, query_title, query_seq.front());
		else
			aligned_targets = TLS_FIX_S390X_MOVE(v);

		if (n == 0 || !new_hits) {
			if (config.query_memory && current_chunk_size >= chunk_size)
				memory->update_failed_count(query_id, current_chunk_size, (i1 - 1)->score);
			const double tail_bit_score = score_matrix.bitscore((i1 - 1)->score);
			if (tail_score == 0 || double((i1 - 1)->score) / (double)tail_score <= config.ranking_score_drop_factor || tail_bit_score < config.ranking_cutoff_bitscore)
				break;
		} else
			tail_score = (i1 - 1)->score;

		i0 = i1;
		i1 = std::min(i1 + std::min(chunk_size, MAX_CHUNK_SIZE), target_scores.cend());
	}

	/*if (multiplier > 1)
		stat.inc(Statistics::HARD_QUERIES);*/

	timer.go("Computing culling");
	culling(aligned_targets, source_query_len, query_title, query_seq.front(), 0);
	if(config.query_memory)
		memory->update(query_id, aligned_targets.begin(), aligned_targets.end());
	stat.inc(Statistics::TARGET_HITS4, aligned_targets.size());
	timer.finish();

	vector<Match> matches = align(aligned_targets, query_seq.data(), query_cb.data(), source_query_len, flags, stat, first_round_traceback);
	std::sort(matches.begin(), matches.end());

	return matches;
}

}