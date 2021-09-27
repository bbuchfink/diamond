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
#include <mutex>
#include "target.h"
#include "../dp/dp.h"
#include "../stats/score_matrix.h"
#include "../data/reference.h"
#include "../util/parallel/thread_pool.h"
#include "../dp/scan_diags.h"

using std::mutex;

namespace Extension {

int gapped_filter(const SeedHit &hit, const LongScoreProfile *query_profile, const Sequence &target, int band, int window, std::function<decltype(DP::ARCH_GENERIC::scan_diags128)> f) {	
	const int slen = (int)target.length();
	const int d = std::max(hit.diag() - band / 2, -(slen - 1)),
		j0 = std::max(hit.j - window, 0),
		j1 = std::min(hit.j + window, slen);
	int scores[128];
	f(query_profile[hit.frame], target, d, j0, j1, scores);
	return DP::diag_alignment(scores, band);
}

bool gapped_filter(FlatArray<SeedHit>::ConstIterator begin, FlatArray<SeedHit>::ConstIterator end, const LongScoreProfile *query_profile, uint32_t target_block_id, Statistics &stat, const Search::Config &params) {
	constexpr int window1 = 100, MIN_STAGE2_QLEN = 100;
		
	const int qlen = (int)query_profile->length();
	const Sequence target = params.target->seqs()[target_block_id];
	const int slen = (int)target.length();
	for (FlatArray<SeedHit>::ConstIterator hit = begin; hit < end; ++hit) {
		stat.inc(Statistics::GAPPED_FILTER_HITS1);
		const int f1 = gapped_filter(*hit, query_profile, target, 64, window1, DP::scan_diags64);
		if (f1 > params.cutoff_gapped1_new(qlen, slen)) {
			stat.inc(Statistics::GAPPED_FILTER_HITS2);
			if (qlen < MIN_STAGE2_QLEN && align_mode.query_translated)
				return true;
			const int f2 = gapped_filter(*hit, query_profile, target, 128, config.gapped_filter_window, DP::scan_diags128);
			if (f2 > params.cutoff_gapped2_new(qlen, slen))
				return true;
		}
	}
	return false;
}

void gapped_filter_worker(size_t i, size_t thread_id, const LongScoreProfile *query_profile, const FlatArray<SeedHit>* seed_hits, const uint32_t* target_block_ids, FlatArray<SeedHit>* out, vector<uint32_t> *target_ids_out, mutex* mtx, const Search::Config *params) {
	thread_local Statistics stat;
	if (gapped_filter(seed_hits->begin(i), seed_hits->end(i), query_profile, target_block_ids[i], stat, *params)) {
		std::lock_guard<mutex> guard(*mtx);
		target_ids_out->push_back(target_block_ids[i]);
		out->push_back(seed_hits->begin(i), seed_hits->end(i));
	}
}

void gapped_filter(const Sequence* query, const Bias_correction* query_cbs, FlatArray<SeedHit>& seed_hits, std::vector<uint32_t>& target_block_ids, Statistics& stat, DP::Flags flags, const Search::Config &params) {
	if (seed_hits.size() == 0)
		return;
	vector<LongScoreProfile> query_profile;
	query_profile.reserve(align_mode.query_contexts);
	for (int i = 0; i < align_mode.query_contexts; ++i)
		if(Stats::CBS::hauser(config.comp_based_stats))
			query_profile.emplace_back(query[i], query_cbs[i]);
		else
			query_profile.emplace_back(query[i]);
	
	FlatArray<SeedHit> hits_out;
	vector<uint32_t> target_ids_out;
	
	if(flag_any(flags, DP::Flags::PARALLEL)) {
		mutex mtx;
		Util::Parallel::scheduled_thread_pool_auto(config.threads_, seed_hits.size(), gapped_filter_worker, query_profile.data(), &seed_hits, target_block_ids.data(), &hits_out, &target_ids_out, &mtx, &params);
	}
	else {

		for (size_t i = 0; i < seed_hits.size(); ++i) {
			if (gapped_filter(seed_hits.begin(i), seed_hits.end(i), query_profile.data(), target_block_ids[i], stat, params)) {
				target_ids_out.push_back(target_block_ids[i]);
				hits_out.push_back(seed_hits.begin(i), seed_hits.end(i));
			}
		}

	}

	seed_hits = std::move(hits_out);
	target_block_ids = std::move(target_ids_out);
}

}