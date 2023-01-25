/****
DIAMOND protein aligner
Copyright (C) 2020-2021 Max Planck Society for the Advancement of Science e.V.

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
#include <numeric>
#include "extend.h"
#include "../data/queries.h"
#include "../basic/config.h"
#include "../stats/hauser_correction.h"
#include "target.h"
#include "../dp/dp.h"
#include "../util/log_stream.h"
#include "../data/reference.h"
#include "../util/system.h"
#include "../util/util.h"
#include "global_ranking/global_ranking.h"
#include "../masking/masking.h"
#include "../search/hit.h"
#include "load_hits.h"
#include "def.h"

using std::accumulate;
using std::vector;
using std::list;
using std::array;
using std::pair;
using std::endl;
using std::move;
using std::make_move_iterator;
using std::any_of;
using std::numeric_limits;
using std::tie;
using std::min;

const SEMap<Extension::Mode> EnumTraits<Extension::Mode>::from_string = {
	{ "banded-fast", Extension::Mode::BANDED_FAST},
	{ "banded-slow", Extension::Mode::BANDED_SLOW},
	{ "full", Extension::Mode::FULL},
	{ "global", Extension::Mode::GLOBAL}
};

namespace Extension {

#include "extend_chunk.h"

const std::map<Sensitivity, Mode> default_ext_mode = {
	{ Sensitivity::FASTER, Mode::BANDED_FAST},
	{ Sensitivity::FAST, Mode::BANDED_FAST},
	{ Sensitivity::DEFAULT, Mode::BANDED_FAST},
	{ Sensitivity::MID_SENSITIVE, Mode::BANDED_FAST},
	{ Sensitivity::SENSITIVE, Mode::BANDED_FAST},
	{ Sensitivity::MORE_SENSITIVE, Mode::BANDED_SLOW},
	{ Sensitivity::VERY_SENSITIVE, Mode::BANDED_SLOW},
	{ Sensitivity::ULTRA_SENSITIVE, Mode::BANDED_SLOW}
};

constexpr int64_t MAX_CHUNK_SIZE = 400, MIN_CHUNK_SIZE = 128, MAPANY_CHUNK_SIZE = 16;

int64_t ranking_chunk_size(int64_t target_count, const int64_t ref_letters, const int64_t max_target_seqs) {
	if (config.no_ranking || config.global_ranking_targets > 0)
		return target_count;
	if (config.ext_chunk_size > 0)
		return config.ext_chunk_size;
	if (config.mapany)
		return MAPANY_CHUNK_SIZE;
	const double default_letters = config.sensitivity >= Sensitivity::VERY_SENSITIVE ? 800 * 1e6 : 2 * 1e9;
	const int64_t block_mult = std::max(int64_t(std::round((double)ref_letters / default_letters)), (int64_t)1);
	if (config.toppercent < 100.0)
		return MIN_CHUNK_SIZE * block_mult;
	const int64_t size = std::max(MIN_CHUNK_SIZE, std::min(make_multiple(max_target_seqs, (int64_t)32), MAX_CHUNK_SIZE)) * block_mult;
	return config.target_hard_cap ? std::min(size, config.target_hard_cap) : size;
}

static bool have_filters(const Search::Config& cfg) {
	return config.min_id > 0 || config.approx_min_id.get(0.0) > 0 || config.query_cover > 0 || config.subject_cover > 0 || config.query_or_target_cover > 0;
}

static HspValues first_round_hspv(const Search::Config& cfg) {
	HspValues first_round = HspValues::NONE;
	if (config.min_id > 0)
		first_round |= HspValues::IDENT | HspValues::LENGTH;
	if (config.query_cover > 0)
		first_round |= HspValues::QUERY_COORDS;
	if (config.subject_cover > 0)
		first_round |= HspValues::TARGET_COORDS;
	if (config.cluster_threshold.present())
		first_round |= cfg.output_format->hsp_values;
	return first_round;
}

static bool ranking_terminate(bool new_hits, int last_tail_score, int tail_score, int64_t targets_processed, int64_t targets_aligned) {
	if (config.target_hard_cap && targets_processed >= config.target_hard_cap)
		return true;
	if (config.mapany && config.toppercent == 100.0 && targets_aligned > 0)
		return true;
	return !new_hits && (last_tail_score == 0
		|| double(tail_score) / (double)last_tail_score <= config.ranking_score_drop_factor
		|| score_matrix.bitscore(tail_score) < config.ranking_cutoff_bitscore);
}

Match Match::self_match(BlockId query_id, Sequence query_seq) {
	Match m(query_id, query_seq, ::Stats::TargetMatrix(), 0, numeric_limits<Score>::max(), 0.0);
	m.hsp.emplace_back();
	m.hsp.back().evalue = 0.0;
	m.hsp.back().score = numeric_limits<Score>::max();
	m.hsp.back().bit_score = DBL_MAX;
	m.hsp.back().query_range = { 0,query_seq.length() };
	m.hsp.back().query_source_range = { 0,query_seq.length() };
	m.hsp.back().subject_range = { 0,query_seq.length() };
	return m;
}

static bool add_self_aln(const Search::Config& cfg) {
	return config.add_self_aln && ((config.self && cfg.current_ref_block == 0) || (!config.self && cfg.current_query_block == cfg.current_ref_block));
}

pair<vector<Match>, Stats> extend(
	BlockId query_id,
	const Search::Config& cfg,
	Statistics& stat,
	DP::Flags flags,
	SeedHitList& l)
{
	const unsigned UNIFIED_TARGET_LEN = 50;
	const unsigned contexts = align_mode.query_contexts;
	vector<Sequence> query_seq;
	vector<Bias_correction> query_cb;
	const char* query_title = cfg.query->ids()[query_id];

	if (config.log_query || (flag_any(flags, DP::Flags::PARALLEL) && !config.swipe_all))
		log_stream << "Query=" << query_title << " Hits=" << l.seed_hits.data_size() << endl;

	for (unsigned i = 0; i < contexts; ++i)
		query_seq.push_back(cfg.query->seqs()[query_id * contexts + i]);
	const unsigned query_len = (unsigned)query_seq.front().length();

	if (::Stats::CBS::hauser(config.comp_based_stats)) {
		for (unsigned i = 0; i < contexts; ++i)
			query_cb.emplace_back(query_seq[i]);
	}
	::Stats::Composition query_comp;
	if (::Stats::CBS::matrix_adjust(config.comp_based_stats))
		query_comp = ::Stats::composition(query_seq[0]);

	const int source_query_len = align_mode.query_translated ? (int)cfg.query->source_seqs()[query_id].length() : (int)cfg.query->seqs()[query_id].length();
	const double self_aln_score = cfg.query->has_self_aln() ? cfg.query->self_aln_score(query_id) : 0.0;
	const size_t target_count = l.target_block_ids.size();
	const int64_t chunk_size = ranking_chunk_size(target_count, cfg.target->seqs().letters(), cfg.max_target_seqs);
	vector<TargetScore>::const_iterator i0 = l.target_scores.cbegin(), i1 = i0 + std::min((ptrdiff_t)chunk_size, l.target_scores.cend() - i0);

	if (config.toppercent == 100.0 && config.min_bit_score == 0.0 && (i1 - i0) < cfg.max_target_seqs && (config.ext_chunk_size == 0 || config.lin_stage1))
#ifdef EVAL_TARGET
		while (i1 < l.target_scores.cend() && i1->evalue <= config.max_evalue && size_t(i1 - i0) < config.max_alignments) ++i1;
#else
		while (i1 < l.target_scores.cend() && score_matrix.evalue(i1->score, query_len, UNIFIED_TARGET_LEN) <= config.max_evalue) i1 += min((ptrdiff_t)16, l.target_scores.cend() - i1);
#endif
	
	const HspValues first_round_hspv = (config.prefix_scan || config.anchored_swipe) ? HspValues::COORDS : HspValues::NONE;
	const bool first_round_culling = !have_filters(cfg) || config.toppercent != 100.0;
	bool new_hits_ev = false;
	int tail_score = 0, previous_tail_score = 0;
	FlatArray<SeedHit> seed_hits_chunk;
	vector<uint32_t> target_block_ids_chunk;
	vector<Match> matches;
	Stats stats;

	do {

		vector<Target> aligned_targets;
		bool new_hits;

		do {
			const int64_t current_chunk_size = i1 - i0;
			const bool multi_chunk = current_chunk_size < (int64_t)l.target_scores.size();

			if (multi_chunk) {
				target_block_ids_chunk.clear();
				seed_hits_chunk.clear();
				target_block_ids_chunk.reserve(i1 - i0);
				seed_hits_chunk.reserve(i1 - i0, accumulate(i0, i1, (int64_t)0, [&l](int64_t n, const TargetScore& s) { return l.seed_hits.count(s.target) + n; }));
				for (vector<TargetScore>::const_iterator j = i0; j < i1; ++j) {
					target_block_ids_chunk.push_back(l.target_block_ids[j->target]);
					seed_hits_chunk.push_back(l.seed_hits.begin(j->target), l.seed_hits.end(j->target));
				}
			}

			pair<vector<Target>, Stats> v = extend(
				query_id,
				query_seq.data(),
				source_query_len,
				query_cb.data(),
				query_comp,
				multi_chunk ? seed_hits_chunk.begin() : l.seed_hits.begin(),
				multi_chunk ? seed_hits_chunk.end() : l.seed_hits.end(),
				multi_chunk ? target_block_ids_chunk.cbegin() : l.target_block_ids.cbegin(),
				cfg,
				stat,
				flags,
				HspValues::NONE);

			stats += v.second;
			stat.inc(Statistics::TARGET_HITS4, v.first.size());
			new_hits = new_hits_ev = v.first.size() > 0;
			if (multi_chunk)
				new_hits = append_hits(aligned_targets, v.first.begin(), v.first.end(), first_round_culling, cfg);
			else
				aligned_targets = move(v.first);

			i0 = i1;
			i1 += std::min(std::min(chunk_size, MAX_CHUNK_SIZE), l.target_scores.cend() - i1);
			previous_tail_score = tail_score;
			if (new_hits)
				tail_score = (i1 - 1)->score;
		} while (i0 < l.target_scores.cend() && !ranking_terminate(new_hits, previous_tail_score, (i1 - 1)->score, i1 - l.target_scores.cbegin(), aligned_targets.size()));

		if (config.swipe_all)
			aligned_targets = full_db_align(query_seq.data(), query_cb.data(), flags, HspValues::NONE, stat, *cfg.target);

		culling(aligned_targets, !first_round_culling, cfg);
		stat.inc(Statistics::TARGET_HITS5, aligned_targets.size());
		
		vector<Match> round_matches = align(aligned_targets, matches.size(), query_seq.data(), query_title, query_cb.data(), source_query_len, self_aln_score, flags, first_round_hspv, first_round_culling, stat, cfg);
		matches.insert(matches.end(), make_move_iterator(round_matches.begin()), make_move_iterator(round_matches.end()));
	} while (config.toppercent == 100.0 && (int64_t)matches.size() < config.max_target_seqs_.get(DEFAULT_MAX_TARGET_SEQS) && i0 < l.target_scores.cend() && new_hits_ev && (!config.mapany || (config.mapany && matches.empty())));

	if (add_self_aln(cfg) && !any_of(matches.cbegin(), matches.cend(), [query_id](const Match& m) { return m.target_block_id == query_id; }))
		matches.push_back(Match::self_match(query_id, query_seq[0]));

	culling(matches, cfg);

	return make_pair(move(matches), stats);
}

pair<vector<Match>, Stats> extend(BlockId query_id, Search::Hit* begin, Search::Hit* end, const Search::Config &cfg, Statistics &stat, DP::Flags flags) {
	task_timer timer(flag_any(flags, DP::Flags::PARALLEL) ? config.target_parallel_verbosity : UINT_MAX);
	timer.go("Loading seed hits");
	SeedHitList l = load_hits(begin, end, cfg.target->seqs());
	stat.inc(Statistics::TARGET_HITS0, l.target_block_ids.size());
	stat.inc(Statistics::TIME_LOAD_HIT_TARGETS, timer.microseconds());
	timer.finish();

	vector<Match> trivial_matches;
	if (config.filter_kmer_len) {
		tie(l, trivial_matches) = kmer_filter(cfg.query->seqs()[query_id], Bias_correction(cfg.query->seqs()[query_id]).int8.data(), *cfg.target, l);
		stat.inc(Statistics::TRIVIAL_ALN, trivial_matches.size());
	}

	const int64_t target_count = (int64_t)l.target_block_ids.size();
	if (target_count == 0 && !config.swipe_all) {
		if (add_self_aln(cfg))
			return { {Match::self_match(query_id, cfg.query->seqs()[query_id])}, Stats() };
		culling(trivial_matches, cfg);
		return { trivial_matches, Stats() };
	}
	const int64_t chunk_size = ranking_chunk_size(target_count, cfg.target->seqs().letters(), cfg.max_target_seqs);

	if (chunk_size < target_count || config.global_ranking_targets > 0) {
		timer.go("Sorting targets by score");
		std::sort(l.target_scores.begin(), l.target_scores.end());
		stat.inc(Statistics::TIME_SORT_TARGETS_BY_SCORE, timer.microseconds());
		timer.finish();
		if (config.global_ranking_targets > 0)
			return make_pair(GlobalRanking::ranking_list(query_id, l.target_scores.begin(), l.target_scores.end(), l.target_block_ids.begin(), l.seed_hits, cfg), Stats());
	}

	pair<vector<Match>, Stats> r = extend(query_id, cfg, stat, flags, l);
	if (!trivial_matches.empty()) {
		r.first.insert(r.first.end(), make_move_iterator(trivial_matches.begin()), make_move_iterator(trivial_matches.end()));
		culling(r.first, cfg);
	}
	return r;
}

}