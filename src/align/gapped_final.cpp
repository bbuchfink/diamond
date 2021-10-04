#include <array>
#include "target.h"
#include "../dp/dp.h"
#include "../output/output_format.h"

using std::array;
using std::list;

namespace Extension {

static void add_dp_targets(const Target& target, int target_idx, const Sequence* query_seq, array<DP::Targets, MAX_CONTEXT>& dp_targets, DP::Flags flags, const HspValues hsp_values, const Mode mode) {
	const Stats::TargetMatrix* matrix = target.adjusted_matrix() ? &target.matrix : nullptr;
	const Loc tlen = target.seq.length();
	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		const Loc qlen = query_seq[frame].length();
		for (const Hsp& hsp : target.hsp[frame]) {
			const int64_t dp_size = flag_any(flags, DP::Flags::FULL_MATRIX)
				? (int64_t)qlen * (int64_t)tlen
				: (int64_t)DpTarget::banded_cols(qlen, tlen, hsp.d_begin, hsp.d_end) * int64_t(hsp.d_end - hsp.d_begin);
			const int b = DP::BandedSwipe::bin(hsp_values, flag_any(flags, DP::Flags::FULL_MATRIX) ? qlen : hsp.d_end - hsp.d_begin, hsp.score, 0, dp_size, matrix ? matrix->score_width() : 0, 0);
			dp_targets[frame][b].emplace_back(target.seq, tlen, hsp.d_begin, hsp.d_end, target_idx, qlen, matrix);
		}
	}
}

vector<Match> align(vector<Target>& targets, const Sequence* query_seq, const Bias_correction* query_cb, int source_query_len, DP::Flags flags, const HspValues first_round, const Mode mode, Statistics& stat) {
	array<DP::Targets, MAX_CONTEXT> dp_targets;
	vector<Match> r;
	if (targets.empty())
		return r;
	r.reserve(targets.size());

	HspValues hsp_values = output_format->hsp_values;
	if (config.max_hsps == 1 && flag_all(first_round, hsp_values)) {
		for (Target& t : targets)
			r.emplace_back(t.block_id, t.seq, t.matrix, t.hsp, t.ungapped_score);
		return r;
	}

	if (mode == Mode::FULL)
		flags |= DP::Flags::FULL_MATRIX;
	if (mode == Mode::GLOBAL)
		flags |= DP::Flags::SEMI_GLOBAL;
	if (config.max_hsps != 1)
		hsp_values |= HspValues::QUERY_COORDS | HspValues::TARGET_COORDS;

	for (size_t i = 0; i < targets.size(); ++i) {
		/*if (config.log_subject)
			std::cout << "Target=" << ref_ids::get()[targets[i].block_id] << " id=" << i << endl;*/
		add_dp_targets(targets[i], i, query_seq, dp_targets, flags, hsp_values, mode);
		r.emplace_back(targets[i].block_id, targets[i].seq, targets[i].matrix, targets[i].ungapped_score);
	}

	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (dp_targets[frame].empty())
			continue;
		DP::Params params{
			query_seq[frame],
			Frame(frame),
			source_query_len,
			Stats::CBS::hauser(config.comp_based_stats) ? query_cb[frame].int8.data() : nullptr,
			flags,
			hsp_values,
			stat
		};
		list<Hsp> hsp = DP::BandedSwipe::swipe(dp_targets[frame], params);
		while (!hsp.empty())
			r[hsp.front().swipe_target].add_hit(hsp, hsp.begin());
	}

	for (Match& match : r)
		match.inner_culling();

	recompute_alt_hsps(r, query_seq, source_query_len, query_cb, hsp_values, stat);

	return r;
}

}
