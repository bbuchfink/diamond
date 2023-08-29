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

#include <vector>
#include "target.h"
#include "../search/hit.h"
#include "../data/sequence_set.h"

namespace Extension {

template<typename It>
static int64_t count_targets(It begin, It end) {
	if (begin == end)
		return 0;
	auto last = (begin++)->subject_;
	int64_t n = 1;
	for (It i = begin; i < end; ++i)
		if (i->subject_ != last) {
			last = i->subject_;
			++n;
		}
	return n;
}

template<typename It>
static SeedHitList load_hits(It begin, It end, const SequenceSet& ref_seqs) {
	std::sort(begin, end, Search::Hit::CmpSubject());	
	const auto targets = count_targets(begin, end);
	const auto hits = end - begin;
	SeedHitList list;
	list.seed_hits.reserve(targets, hits);
	list.target_block_ids.reserve(targets);
	list.target_scores.reserve(targets);

	if (hits <= 0)
		return list;
	
	const size_t total_subjects = ref_seqs.size();
	unsigned target_len;
	uint32_t target = UINT32_MAX;
	uint16_t score = 0;
#ifdef HIT_KEEP_TARGET_ID
	if(true) {
#else
	if (std::log2(total_subjects) * hits < total_subjects / 10) {
#endif
		for (auto i = begin; i < end; ++i) {
#ifdef HIT_KEEP_TARGET_ID
			std::pair<size_t, size_t> l{ (size_t)i->target_block_id, (size_t)i->subject_ - ref_seqs.position(i->target_block_id, 0) };
#else
			std::pair<size_t, size_t> l = ref_seqs.local_position((uint64_t)i->subject_);
#endif
			const uint32_t t = (uint32_t)l.first;
			if (t != target) {
				if (target != UINT32_MAX) {
#ifdef EVAL_TARGET
					list.target_scores.push_back({ uint32_t(list.target_block_ids.size() - 1), score, score_matrix.evalue(score, query_len, target_len) });
#else
					list.target_scores.push_back({ uint32_t(list.target_block_ids.size() - 1), score });
#endif
					score = 0;
				}
				list.seed_hits.next();
				target = t;
				target_len = (unsigned)ref_seqs[target].length();
				list.target_block_ids.push_back(target);
			}
			list.seed_hits.push_back({ (int)i->seed_offset_, (int)l.second, i->score_, (unsigned)i->query_ % align_mode.query_contexts });
			score = std::max(score, i->score_);
		}
	}
	else {
		typename std::vector<int64_t>::const_iterator limit_begin = ref_seqs.limits_begin(), it = limit_begin;
		for (auto i = begin; i < end; ++i) {
			const int64_t subject_offset = i->subject_;
			while (*it <= subject_offset) ++it;
			uint32_t t = (uint32_t)(it - limit_begin) - 1;
			if (t != target) {
				if (target != UINT32_MAX) {
#ifdef EVAL_TARGET
					list.target_scores.push_back({ uint32_t(list.target_block_ids.size() - 1), score, score_matrix.evalue(score, query_len, target_len) });
#else
					list.target_scores.push_back({ uint32_t(list.target_block_ids.size() - 1), score });
#endif
					score = 0;
				}
				list.seed_hits.next();
				list.target_block_ids.push_back(t);
				target = t;
				target_len = (unsigned)ref_seqs[target].length();
			}
			list.seed_hits.push_back({ (int)i->seed_offset_, (int)(subject_offset - *(it - 1)), i->score_, (unsigned)i->query_ % align_mode.query_contexts });
			score = std::max(score, i->score_);
		}
	}
	if (target != UINT32_MAX)
#ifdef EVAL_TARGET
		list.target_scores.push_back({ uint32_t(list.target_block_ids.size() - 1), score, score_matrix.evalue(score, query_len, target_len) });
#else
		list.target_scores.push_back({ uint32_t(list.target_block_ids.size() - 1), score });
#endif

	return list;
}

}