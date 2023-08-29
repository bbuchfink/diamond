/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

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

#include "global_ranking.h"
#define _REENTRANT
#include "../../lib/ips4o/ips4o.hpp"
#include "../../search/hit.h"
#include "../../util/data_structures/deque.h"
#include "../../util/util.h"
#include "../../util/algo/algo.h"
#include "../../data/seed_array.h"
#if _MSC_FULL_VER == 191627042
#include "../../util/algo/merge_sort.h"
#endif
#include "../load_hits.h"
#include "../dp/ungapped.h"

using std::endl;
using std::thread;
using std::pair;
using std::vector;
using SeedHits = Search::Config::RankingBuffer;

// #define BATCH_BINSEARCH

namespace Extension { namespace GlobalRanking {

static void get_query_hits(SeedHits::Iterator begin, SeedHits::Iterator end, vector<Hit>& hits, Search::Config& cfg) {
	hits.clear();
	const SequenceSet& target_seqs = cfg.target->seqs();
#ifdef KEEP_TARGET_ID
	auto get_target = [](const Search::Hit& hit) { return (BlockId)hit.subject_; };
	auto it = merge_keys(begin, end, get_target);
	while (it.good()) {
		uint16_t score = 0;
		for (SeedHits::Iterator i = it.begin(); i != it.end(); ++i)
			score = std::max(score, i->score_);
		hits.emplace_back((uint32_t)cfg.target->block_id2oid(it.key()), score, 0);
		++it;
	}
#else
#ifdef BATCH_BINSEARCH
	vector<Hit> hit1;
	hit1.reserve(end - begin);
	target_seqs.local_position_batch(begin, end, std::back_inserter(hit1), Search::Hit::CmpTargetOffset());
	for (size_t i = 0; i < hit1.size(); ++i) {
		hit1[i].score = begin[i].score_;
	}
	auto it = merge_keys(hit1.begin(), hit1.end(), Hit::Target());
	while (it.good()) {
		uint16_t score = 0;
		for (auto i = it.begin(); i != it.end(); ++i)
			score = std::max(score, i->score);
		hits.emplace_back((uint32_t)cfg.target->block_id2oid(it.key()), score);
		++it;
	}
#else
	auto get_target = [&target_seqs](const Search::Hit& hit) { return target_seqs.local_position((uint64_t)hit.subject_).first; };
	auto it = merge_keys(begin, end, get_target);
	while (it.good()) {
		uint16_t score = 0;
		for (SeedHits::Iterator i = it.begin(); i != it.end(); ++i)
			score = std::max(score, i->score_);
		hits.emplace_back((uint32_t)cfg.target->block_id2oid(it.key()), score, 0);
		++it;
	}
#endif
#endif
}

static pair<int, unsigned> target_score(const FlatArray<Extension::SeedHit>::DataIterator begin, const FlatArray<Extension::SeedHit>::DataIterator end, const Sequence* query_seq, const Sequence& target_seq) {
	if (config.no_reextend) {
		int score = begin->score;
		unsigned context = begin->frame;
		for (auto i = begin + 1; i != end; ++i) {
			if (i->score > score) {
				score = i->score;
				context = i->frame;
			}
		}
		return { score,context };
	}
	std::sort(begin, end);
	DiagonalSegment d = xdrop_ungapped(query_seq[begin->frame], nullptr, target_seq, begin->i, begin->j, false);
	int score = d.score;
	unsigned context = begin->frame;
	for (auto i = begin + 1; i != end; ++i) {
		if (d.diag() == i->diag() && d.subject_end() >= i->j)
			continue;
		d = xdrop_ungapped(query_seq[i->frame], nullptr, target_seq, i->i, i->j, false);
		if (d.score > score) {
			score = d.score;
			context = i->frame;
		}
	}
	return { score,context };
}

static void get_query_hits_reextend(size_t source_query_block_id, SeedHits::Iterator begin, SeedHits::Iterator end, vector<Hit>& hits, Search::Config& cfg) {
	const unsigned contexts = align_mode.query_contexts;
	vector<Sequence> query_seq;
	for (unsigned i = 0; i < contexts; ++i)
		query_seq.push_back(cfg.query->seqs()[source_query_block_id * contexts + i]);
	const SequenceSet& target_seqs = cfg.target->seqs();

	hits.clear();
	Extension::SeedHitList l = Extension::load_hits(begin, end, cfg.target->seqs());
	for (size_t i = 0; i < l.target_block_ids.size(); ++i) {
		const auto score = target_score(l.seed_hits.begin(i), l.seed_hits.end(i), query_seq.data(), target_seqs[l.target_block_ids[i]]);
		hits.emplace_back((uint32_t)cfg.target->block_id2oid(l.target_block_ids[i]), score.first, score.second);
	}
}

static void merge_hits(const size_t query, vector<Hit>& hits, vector<Hit>& merged, Search::Config& cfg, size_t& merged_count) {
	const size_t N = config.global_ranking_targets;
	//std::sort(hits.begin(), hits.end());
	vector<Hit>::iterator table_begin = cfg.ranking_table->begin() + query * N, table_end = table_begin + N;
	while (table_end > table_begin && (table_end - 1)->score == 0) --table_end;
	const auto count = table_end - table_begin;
	hits.insert(hits.end(), table_begin, table_end);
	std::sort(hits.begin(), hits.end(), Hit::CmpOidScore());
	merged.clear();
	std::unique_copy(hits.begin(), hits.end(), std::back_inserter(merged), Hit::CmpOid());
	std::sort(merged.begin(), merged.end());
	const auto n = std::min(N, merged.size());
	std::copy(merged.begin(), merged.begin() + n, table_begin);
	merged_count += n - count;
	//merged.clear();
	//merged_count += Util::Algo::merge_capped(table_begin, table_end, hits.begin(), hits.end(), N, std::back_inserter(merged));
	//std::copy(merged.begin(), merged.end(), table_begin);
}

void update_table(Search::Config& cfg) {
	SeedHits& hits = *cfg.global_ranking_buffer;
	log_stream << "Seed hits = " << hits.size() << endl;
	if (hits.size() == 0)
		return;
	TaskTimer timer("Sorting seed hits");
#if _MSC_FULL_VER == 191627042
	merge_sort(hits.begin(), hits.end(), config.threads_, Search::Hit::CmpQueryTarget());
#else
	ips4o::parallel::sort(hits.begin(), hits.end(), Search::Hit::CmpQueryTarget(), config.threads_);
#endif
	timer.go("Creating partition");
	auto p = Util::Algo::partition_table(hits.begin(), hits.end(), config.threads_ * 8, ::Search::Hit::SourceQuery{ align_mode.query_contexts });
	timer.go("Processing seed hits");
	std::atomic_size_t merged_count(0), next(0);
	auto worker = [&cfg, &merged_count, &next, &p]() {
		vector<Hit> hits, merged;
		while (true) {
			auto i = next++;
			if (i + 1 >= p.size())
				break;
			auto begin = p[i], end = p[i + 1];
			auto it = merge_keys(begin, end, ::Search::Hit::SourceQuery{ align_mode.query_contexts });
			size_t n = 0;
			while (it.good()) {
				const size_t query = it.begin()->query_ / align_mode.query_contexts;
				get_query_hits_reextend(query, it.begin(), it.end(), hits, cfg);
				merge_hits(query, hits, merged, cfg, n);
				++it;
			}
			merged_count += n;
		}
	};
	vector<thread> threads;
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(worker);
	for (thread& t : threads)
		t.join();
	timer.go("Deallocating seed hit list");
	cfg.global_ranking_buffer.reset();
	timer.finish();
	log_stream << "Merged targets = " << merged_count << endl;
}

}}