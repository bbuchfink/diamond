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

#include <algorithm>
#include <utility>
#include <vector>
#include "target.h"
#include "search/hit.h"
#include "data/sequence_set.h"

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

static bool hamming_filter_without_boundaries(const SeedHit& hit, const Sequence& query, const Sequence& target, unsigned hamming_filter_id) {
	constexpr int WINDOW_LEFT = 16, WINDOW_RIGHT = 32;
	const int begin = std::max(std::max(-WINDOW_LEFT, -hit.i), -hit.j);
	const int end = std::min(std::min(WINDOW_RIGHT, (int)query.length() - hit.i), (int)target.length() - hit.j);
	const int len = end - begin;
	if (len <= 0)
		return false;
	unsigned id = 0;
	for (int i = begin; i < end; ++i)
		if (query[hit.i + i] == target[hit.j + i])
			++id;
	return id >= hamming_filter_id;
}

static void filter_hamming_boundary_crossings(SeedHitList& list, const Sequence* query_seq, Loc query_len, const SequenceSet& ref_seqs, unsigned hamming_filter_id) {
	if (hamming_filter_id == 0 || list.seed_hits.data_size() == 0)
		return;
#ifndef EVAL_TARGET
	(void)query_len;
#endif

	SeedHitList filtered;
	filtered.seed_hits.reserve(list.target_block_ids.size(), list.seed_hits.data_size());
	filtered.target_block_ids.reserve(list.target_block_ids.size());
	filtered.target_scores.reserve(list.target_scores.size());
	for (size_t i = 0; i < list.target_block_ids.size(); ++i) {
		const BlockId target_block_id = list.target_block_ids[i];
		const Sequence target = ref_seqs[target_block_id];
		bool target_started = false;
		uint16_t score = 0;
		for (auto hit = list.seed_hits.begin(i); hit != list.seed_hits.end(i); ++hit) {
			if (!hamming_filter_without_boundaries(*hit, query_seq[hit->frame], target, hamming_filter_id))
				continue;
			if (!target_started) {
				filtered.seed_hits.next();
				filtered.target_block_ids.push_back(target_block_id);
				target_started = true;
			}
			filtered.seed_hits.push_back(*hit);
			score = std::max(score, (uint16_t)hit->score);
		}
		if (target_started) {
#ifdef EVAL_TARGET
			filtered.target_scores.push_back({ uint32_t(filtered.target_block_ids.size() - 1), score, score_matrix.evalue(score, query_len, target.length()) });
#else
			filtered.target_scores.push_back({ uint32_t(filtered.target_block_ids.size() - 1), score });
#endif
		}
	}
	list = std::move(filtered);
}

}