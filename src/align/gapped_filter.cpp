#include <algorithm>
#include <utility>
#include "target.h"
#include "../dp/dp.h"
#include "../basic/score_matrix.h"

namespace Extension {

static bool score_cutoff(int score, int query_len) {
	if (config.gapped_filter_evalue > 0.0)
		return score_matrix.evalue(score, query_len) <= config.gapped_filter_evalue;
	else
		return score_matrix.bitscore(score) >= config.gapped_filter_score;
}

static int diag_alignment(const int* s) {
	int best = 0, best_gap = -score_matrix.gap_open(), d = -1;
	for (int i = 0; i < DP::GAPPED_FILTER_BAND; ++i) {
		if (s[i] < config.gapped_filter_diag_score)
			continue;
		const int gap_score = -score_matrix.gap_extend() * (i - d) + best_gap;
		int n = s[i];
		if (gap_score + s[i] > best) {
			best = n = gap_score + s[i];
		}
		if (s[i] > best) {
			best = n = s[i];
		}
		const int open_score = -score_matrix.gap_open() + n;
		if (open_score > gap_score) {
			best_gap = open_score;
			d = i;
		}
	}
	return best;
}

bool gapped_filter(const LongScoreProfile *query_profile, const WorkTarget& target) {
	const int slen = (int)target.seq.length(), qlen = (int)query_profile[0].length();
	int scores[DP::GAPPED_FILTER_BAND];
	for (int frame = 0; frame < align_mode.query_contexts; ++frame)
		for (const Hsp_traits& hsp : target.hsp[frame]) {
			const int d = std::max((hsp.d_max + hsp.d_min) / 2 - DP::GAPPED_FILTER_BAND / 2, -(slen - 1)),
				j0 = std::max(hsp.subject_range.begin_ - config.gapped_filter_window, 0),
				j1 = std::min(hsp.subject_range.end_ + config.gapped_filter_window, slen);
			DP::scan_diags(query_profile[frame], target.seq, d, j0, j1, scores);
			const int score = diag_alignment(scores);
			if (score_cutoff(score, qlen))
				return true;
		}
	return false;
}

vector<WorkTarget> gapped_filter(const sequence* query, const Bias_correction* query_cbs, std::vector<WorkTarget>& targets) {
	vector<WorkTarget> out;
	vector<LongScoreProfile> query_profile;
	query_profile.reserve(align_mode.query_contexts);
	for (unsigned i = 0; i < align_mode.query_contexts; ++i)
		query_profile.emplace_back(query[i], query_cbs[i]);
	const int qlen = (int)query[0].length();

	for (WorkTarget& target : targets) {
		if (score_cutoff(target.filter_score, qlen) || gapped_filter(query_profile.data(), target))
			out.push_back(std::move(target));
	}

	return out;
}

}