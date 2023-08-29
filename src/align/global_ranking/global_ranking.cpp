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

#include <mutex>
#include "global_ranking.h"
#include "../data/queries.h"
#include "../util/sequence/sequence.h"
#include "../../dp/ungapped.h"
#include "../../util/util.h"

using std::mutex;
using std::vector;

namespace Extension { namespace GlobalRanking {

uint16_t recompute_overflow_scores(FlatArray<SeedHit>::DataConstIterator begin, FlatArray<SeedHit>::DataConstIterator end, size_t query_id, uint32_t target_id, const Search::Config& cfg) {
	const auto query = cfg.query->seqs()[query_id];
	const auto target = cfg.target->seqs()[target_id];
	int score = 0;
	for (auto it = begin; it < end; ++it) {
		if (it->score != UCHAR_MAX)
			continue;
		const Sequence query_clipped = Util::Seq::clip(query.data() + it->i - config.ungapped_window, config.ungapped_window * 2, config.ungapped_window);
		const ptrdiff_t window_left = query.data() + it->i - query_clipped.data();
		const int s = ungapped_window(query_clipped.data(), target.data() + it->j - window_left, (int)query_clipped.length());
		score = std::max(score, s);
	}
	return (uint16_t)std::min(score, USHRT_MAX);
}

std::vector<Extension::Match> ranking_list(size_t query_id, std::vector<TargetScore>::iterator begin, std::vector<TargetScore>::iterator end, std::vector<uint32_t>::const_iterator target_block_ids, const FlatArray<SeedHit>& seed_hits, const Search::Config& cfg) {
	size_t overflows = 0;
	for (auto it = begin; it < end && it->score >= UCHAR_MAX; ++it)
		if (it->score == UCHAR_MAX) {
			it->score = recompute_overflow_scores(seed_hits.begin(it->target), seed_hits.end(it->target), query_id, target_block_ids[it->target], cfg);
			++overflows;
		}
	if (overflows > 0)
		std::sort(begin, end); // should also sort by target block id

	int64_t n = 0;
	vector<Match> r;
	r.reserve(std::min(end - begin, (ptrdiff_t)config.global_ranking_targets));
	for (auto i = begin; i < end && n < config.global_ranking_targets; ++i, ++n) {
		r.emplace_back(target_block_ids[i->target], Sequence(), ::Stats::TargetMatrix(), i->score);
	}
	return r;
}

size_t write_merged_query_list_intro(uint32_t query_id, TextBuffer& buf) {
	size_t seek_pos = buf.size();
	buf.write(query_id).write((uint32_t)0);
	return seek_pos;
}

void write_merged_query_list(const IntermediateRecord& r, TextBuffer& out, BitVector& ranking_db_filter, Statistics& stat) {
	/*out.write(r.subject_oid);
	out.write(uint16_t(r.score));
	ranking_db_filter.set(r.subject_oid);
	stat.inc(Statistics::TARGET_HITS1);*/
}

void finish_merged_query_list(TextBuffer& buf, size_t seek_pos) {
	*(uint32_t*)(&buf[seek_pos + sizeof(uint32_t)]) = safe_cast<uint32_t>(buf.size() - seek_pos - sizeof(uint32_t) * 2);
}

QueryList fetch_query_targets(InputFile& query_list, uint32_t& next_query) {
	static mutex mtx;
	std::lock_guard<mutex> lock(mtx);
	QueryList r;
	r.last_query_block_id = next_query;
	uint32_t size;
	try {
		query_list >> r.query_block_id;
	}
	catch (EndOfStream&) {
		return r;
	}
	next_query = r.query_block_id + 1;
	query_list >> size;
	size_t n = size / 6;
	r.targets.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		uint32_t target;
		uint16_t score;
		query_list >> target >> score;
		r.targets.push_back({ target,score });
	}
	return r;
}

}}