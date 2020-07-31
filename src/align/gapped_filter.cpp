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
#include "../basic/score_matrix.h"
#include "../data/reference.h"
#include "../util/parallel/thread_pool.h"
#include "../dp/scan_diags.h"

using std::mutex;

namespace Extension {

/*bool gapped_filter(const LongScoreProfile *query_profile, const WorkTarget& target, Statistics &stat) {
	const int slen = (int)target.seq.length(), qlen = (int)query_profile[0].length();
	int scores[128];
	stat.inc(Statistics::GAPPED_FILTER_TARGETS);
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame)
		for (const Hsp_traits& hsp : target.hsp[frame]) {
			stat.inc(Statistics::GAPPED_FILTER_HITS1);
			const int d = std::max((hsp.d_max + hsp.d_min) / 2 - 128 / 2, -(slen - 1)),
				j0 = std::max(hsp.subject_range.begin_ - config.gapped_filter_window, 0),
				j1 = std::min(hsp.subject_range.end_ + config.gapped_filter_window, slen);
			DP::scan_diags128(query_profile[frame], target.seq, d, j0, j1, scores);
			const int score = DP::diag_alignment(scores, 128);
			if (score_cutoff(score, qlen))
				return true;
		}
	return false;
}

vector<WorkTarget> gapped_filter(const sequence* query, const Bias_correction* query_cbs, std::vector<WorkTarget>& targets, Statistics &stat) {
	vector<WorkTarget> out;
	vector<LongScoreProfile> query_profile;
	query_profile.reserve(align_mode.query_contexts);
	for (unsigned i = 0; i < align_mode.query_contexts; ++i)
		query_profile.emplace_back(query[i], query_cbs[i]);
	const int qlen = (int)query[0].length();

	for (WorkTarget& target : targets) {
		if (score_cutoff(target.filter_score, qlen)) {
			out.push_back(std::move(target));
			continue;
		}
		if(gapped_filter(query_profile.data(), target, stat))
			out.push_back(std::move(target));
	}

	return out;
}*/

int gapped_filter(const SeedHit &hit, const LongScoreProfile *query_profile, const sequence &target, int band, int window, std::function<decltype(DP::ARCH_GENERIC::scan_diags128)> f) {	
	const int slen = (int)target.length();
	const int d = std::max(hit.diag() - band / 2, -(slen - 1)),
		j0 = std::max(hit.j - window, 0),
		j1 = std::min(hit.j + window, slen);
	int scores[128];
	f(query_profile[hit.frame], target, d, j0, j1, scores);
	return DP::diag_alignment(scores, band);
}

bool gapped_filter(const SeedHit *begin, const SeedHit *end, const LongScoreProfile *query_profile, uint32_t target_block_id, Statistics &stat, const Parameters &params) {
	constexpr int window1 = 100;
		
	const int qlen = (int)query_profile->length();
	const sequence target = ref_seqs::get()[target_block_id];
	for (const SeedHit* hit = begin; hit < end; ++hit) {
		stat.inc(Statistics::GAPPED_FILTER_HITS1);
		const int f1 = gapped_filter(*hit, query_profile, target, 64, window1, DP::scan_diags64);
		if(f1 > params.cutoff_gapped1(qlen)) {
			stat.inc(Statistics::GAPPED_FILTER_HITS2);
			const int f2 = gapped_filter(*hit, query_profile, target, 128, config.gapped_filter_window, DP::scan_diags128);
			if(f2 > params.cutoff_gapped2(qlen))
				return true;
		}
	}
	return false;
}

void gapped_filter_worker(size_t i, size_t thread_id, const LongScoreProfile *query_profile, const FlatArray<SeedHit>* seed_hits, const uint32_t* target_block_ids, FlatArray<SeedHit>* out, vector<uint32_t> *target_ids_out, mutex* mtx, const Parameters *params) {
	thread_local Statistics stat;
	if (gapped_filter(seed_hits->begin(i), seed_hits->end(i), query_profile, target_block_ids[i], stat, *params)) {
		std::lock_guard<mutex> guard(*mtx);
		target_ids_out->push_back(target_block_ids[i]);
		out->push_back(seed_hits->begin(i), seed_hits->end(i));
	}
}

void gapped_filter(const sequence* query, const Bias_correction* query_cbs, FlatArray<SeedHit>& seed_hits, std::vector<uint32_t>& target_block_ids, Statistics& stat, int flags, const Parameters &params) {
	if (seed_hits.size() == 0)
		return;
	vector<LongScoreProfile> query_profile;
	query_profile.reserve(align_mode.query_contexts);
	for (unsigned i = 0; i < align_mode.query_contexts; ++i)
		if(config.comp_based_stats)
			query_profile.emplace_back(query[i], query_cbs[i]);
		else
			query_profile.emplace_back(query[i]);
	
	FlatArray<SeedHit> hits_out;
	vector<uint32_t> target_ids_out;
	
	if(flags & TARGET_PARALLEL) {
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