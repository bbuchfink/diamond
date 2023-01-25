#include <map>
#include "target.h"
#include "../dp/dp.h"

using std::list;
using std::map;
using std::vector;

namespace Extension {

vector<Target> full_db_align(const Sequence* query_seq, const Bias_correction* query_cb, DP::Flags flags, const HspValues hsp_values, Statistics& stat, const Block& target_block) {
	vector<DpTarget> v;
	vector<Target> r;
	::Stats::TargetMatrix matrix;
	list<Hsp> hsp;
	const SequenceSet& ref_seqs = target_block.seqs();

	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		DP::Params params{
			query_seq[frame],
			"",
			Frame(frame),
			query_seq[frame].length(),
			::Stats::CBS::hauser(config.comp_based_stats) ? query_cb[frame].int8.data() : nullptr,
			flags | DP::Flags::FULL_MATRIX,
			hsp_values,
			stat,
			nullptr
		};
		list<Hsp> frame_hsp = DP::BandedSwipe::swipe_set(ref_seqs.cbegin(), ref_seqs.cend(), params);
		hsp.splice(hsp.begin(), frame_hsp, frame_hsp.begin(), frame_hsp.end());
	}

	map<BlockId, BlockId> subject_idx;
	while (!hsp.empty()) {
		BlockId block_id = hsp.begin()->swipe_target;
		const auto it = subject_idx.emplace(block_id, (BlockId)r.size());
		if (it.second)
			r.emplace_back(block_id, ref_seqs[block_id], 0, matrix);
		BlockId i = it.first->second;
		r[i].add_hit(hsp, hsp.begin());
	}

	return r;
}

}