#include <array>
#include "target.h"
#include "../dp/dp.h"
#include "../output/output_format.h"
#include "../util/util.h"
#include "def.h"

using std::array;
using std::list;
using std::vector;

namespace Extension {

HspValues filter_hspvalues() {
	HspValues hsp_values = HspValues::NONE;
	if (config.max_hsps != 1)
		hsp_values |= HspValues::QUERY_COORDS | HspValues::TARGET_COORDS;
	if (config.min_id > 0)
		hsp_values |= HspValues::IDENT | HspValues::LENGTH;
	if (config.approx_min_id.get(0.0) > 0)
		hsp_values |= HspValues::COORDS;
	if (config.query_cover > 0)
		hsp_values |= HspValues::QUERY_COORDS;
	if (config.subject_cover > 0)
		hsp_values |= HspValues::TARGET_COORDS;
	if (config.query_or_target_cover > 0)
		hsp_values |= HspValues::COORDS;
	return hsp_values;
}

static bool first_round_filter_all(const Search::Config& cfg, HspValues first_round_hsp_values) {
	if (config.min_id > 0 && !flag_all(first_round_hsp_values, HspValues::IDENT | HspValues::LENGTH))
		return false;
	if (config.approx_min_id.get(0.0) > 0 && !flag_all(first_round_hsp_values, HspValues::COORDS))
		return false;
	if (config.query_cover > 0 && !flag_all(first_round_hsp_values, HspValues::QUERY_COORDS))
		return false;
	if (config.subject_cover > 0 && !flag_all(first_round_hsp_values, HspValues::TARGET_COORDS))
		return false;
	if (config.query_or_target_cover > 0 && !flag_all(first_round_hsp_values, HspValues::COORDS))
		return false;
	return true;
}

static void add_dp_targets(const Target& target, int target_idx, const Sequence* query_seq, array<DP::Targets, MAX_CONTEXT>& dp_targets, DP::Flags flags, const HspValues hsp_values, const Mode mode) {
	const ::Stats::TargetMatrix* matrix = target.adjusted_matrix() ? &target.matrix : nullptr;
	const Loc tlen = target.seq.length();
	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		const Loc qlen = query_seq[frame].length();
		for (const Hsp& hsp : target.hsp[frame]) {
			const int64_t dp_size = flag_any(flags, DP::Flags::FULL_MATRIX)
				? (int64_t)qlen * (int64_t)tlen
				: (int64_t)DpTarget::banded_cols(qlen, tlen, hsp.d_begin, hsp.d_end) * int64_t(hsp.d_end - hsp.d_begin);
			const int b = DP::BandedSwipe::bin(hsp_values, flag_any(flags, DP::Flags::FULL_MATRIX) ? qlen : hsp.d_end - hsp.d_begin, hsp.score, 0, dp_size, matrix ? matrix->score_width() : 0, 0);
			dp_targets[frame][b].emplace_back(target.seq, tlen, hsp.d_begin, hsp.d_end, Interval(), 0, target_idx, qlen, matrix);
		}
	}
}

vector<Match> align(vector<Target>& targets, const int64_t previous_matches, const Sequence* query_seq, const char* query_id, const Bias_correction* query_cb, int source_query_len, double query_self_aln_score, DP::Flags flags, const HspValues first_round, const bool first_round_culling, Statistics& stat, const Search::Config& cfg) {
	static const int64_t MIN_STEP = 16;
	vector<Match> r;
	if (targets.empty())
		return r;
	HspValues hsp_values = cfg.output_format->hsp_values;
	const bool copy_all = config.max_hsps == 1 && flag_all(first_round, hsp_values) && first_round_filter_all(cfg, first_round);
	if (copy_all)
		r.reserve(targets.size());
	for (Target& t : targets)
		if (copy_all || t.done)
			r.emplace_back(t.block_id, t.seq, t.matrix, t.hsp, t.ungapped_score);
	if (r.size() == targets.size()) {
		apply_filters(r.begin(), r.end(), source_query_len, query_id, query_self_aln_score, query_seq[0], cfg);
		return r;
	}
	
	if (cfg.extension_mode == Mode::FULL)
		flags |= DP::Flags::FULL_MATRIX;
	if (cfg.extension_mode == Mode::GLOBAL)
		flags |= DP::Flags::SEMI_GLOBAL;
	hsp_values |= filter_hspvalues();

	vector<Target>::iterator it = targets.begin();
	auto goon = [&r, &cfg, previous_matches]() { return config.toppercent == 100.0 ? ((int64_t)r.size() + previous_matches) < cfg.max_target_seqs : true; };

	do {
		array<DP::Targets, MAX_CONTEXT> dp_targets;
		const int64_t step_size = !first_round_culling && config.toppercent == 100.0
			? std::min(make_multiple(std::max(cfg.max_target_seqs - (int64_t)r.size(), MIN_STEP), MIN_STEP), targets.end() - it)
			: targets.size();
		
		r.reserve(r.size() + step_size);

		const int64_t matches_begin = r.size();

		for (auto i = it; i < it + step_size; ++i) {
			/*if (config.log_subject)
				std::cout << "Target=" << ref_ids::get()[targets[i].block_id] << " id=" << i << endl;*/
			if (i->done)
				continue;
			add_dp_targets(*i, (int32_t)r.size(), query_seq, dp_targets, flags, hsp_values, cfg.extension_mode);
			r.emplace_back(i->block_id, i->seq, i->matrix, i->ungapped_score);
		}

		for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
			if (dp_targets[frame].empty())
				continue;
			DP::Params params{
				query_seq[frame],
				query_id,
				Frame(frame),
				source_query_len,
				::Stats::CBS::hauser(config.comp_based_stats) ? query_cb[frame].int8.data() : nullptr,
				flags,
				hsp_values,
				stat,
				cfg.thread_pool.get()
			};
			list<Hsp> hsp = DP::BandedSwipe::swipe(dp_targets[frame], params);
			while (!hsp.empty())
				r[hsp.front().swipe_target].add_hit(hsp, hsp.begin());
		}

		for (int64_t i = matches_begin; i < (int64_t)r.size(); ++i)
			r[i].inner_culling();

		apply_filters(r.begin() + matches_begin, r.end(), source_query_len, query_id, query_self_aln_score, query_seq[0], cfg);
		culling(r, cfg);

		stat.inc(Statistics::TARGET_HITS6, step_size);
		it += step_size;

	} while (it < targets.end() && goon());

	recompute_alt_hsps(r.begin(), r.end(), query_seq, source_query_len, query_cb, hsp_values, stat);
	return r;
}

}
