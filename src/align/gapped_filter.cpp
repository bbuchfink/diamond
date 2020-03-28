#include <algorithm>
#include <utility>
#include "target.h"
#include "../dp/dp.h"
#include "../basic/score_matrix.h"
#include "../data/reference.h"

namespace Extension {

static bool score_cutoff(int score, int query_len) {
	if (config.gapped_filter_evalue > 0.0) {
		return score_matrix.evalue(score, query_len) <= config.gapped_filter_evalue;
	}
	else
		return score_matrix.bitscore(score) >= config.gapped_filter_score;
}

bool gapped_filter(const LongScoreProfile *query_profile, const WorkTarget& target, Statistics &stat) {
	const int slen = (int)target.seq.length(), qlen = (int)query_profile[0].length();
	int scores[128];
	stat.inc(Statistics::GAPPED_FILTER_TARGETS);
	for (int frame = 0; frame < align_mode.query_contexts; ++frame)
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
}

int gapped_filter(const SeedHit &hit, const LongScoreProfile *query_profile, size_t target_block_id, int band, int window, std::function<decltype(DP::ARCH_GENERIC::scan_diags128)> f) {
	const sequence target = ref_seqs::get()[target_block_id];
	const int slen = (int)target.length();
	const int d = std::max(hit.diag() - band / 2, -(slen - 1)),
		j0 = std::max(hit.j - window, 0),
		j1 = std::min(hit.j + window, slen);
	int scores[128];
	f(query_profile[hit.frame], target, d, j0, j1, scores);
	return DP::diag_alignment(scores, band);
}

void gapped_filter(const sequence* query, const Bias_correction* query_cbs, FlatArray<SeedHit>& seed_hits, std::vector<size_t>& target_block_ids, Statistics& stat) {
	constexpr int window1 = 100;
	constexpr double evalue1 = 1.0e+04;

	if (seed_hits.size() == 0)
		return;
	vector<LongScoreProfile> query_profile;
	query_profile.reserve(align_mode.query_contexts);
	for (unsigned i = 0; i < align_mode.query_contexts; ++i)
		query_profile.emplace_back(query[i], query_cbs[i]);
	const int qlen = (int)query[0].length();

	FlatArray<SeedHit> hits_out;
	vector<size_t> target_ids_out;

	for (size_t i = 0; i < seed_hits.size(); ++i) {
		int n = 0;
		hits_out.next();
		for (const SeedHit* hit = seed_hits.begin(i); hit < seed_hits.end(i); ++hit) {
			stat.inc(Statistics::GAPPED_FILTER_HITS1);
			const int f1 = gapped_filter(*hit, query_profile.data(), target_block_ids[i], 64, window1, DP::scan_diags64);
			if(score_matrix.evalue(f1, qlen) <= evalue1) {
				stat.inc(Statistics::GAPPED_FILTER_HITS2);
				const int f2 = gapped_filter(*hit, query_profile.data(), target_block_ids[i], 128, config.gapped_filter_window, DP::scan_diags128);
				if (score_cutoff(f2, qlen)) {
					hits_out.push_back(*hit);
					++n;
				}
			}
		}
		if (n > 0)
			target_ids_out.push_back(target_block_ids[i]);
		else
			hits_out.pop_back();
	}

	seed_hits = std::move(hits_out);
	target_block_ids = std::move(target_ids_out);
}

}