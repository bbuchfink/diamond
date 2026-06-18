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

#include <map>
#include "target.h"
#include "dp/dp.h"

using std::list;
using std::map;
using std::vector;

namespace Extension {

vector<Target> full_db_align(const Sequence* query_seq, const HauserCorrection* query_cb, DP::Flags flags, const HspValues hsp_values, Statistics& stat, const Block& target_block) {
	vector<DpTarget> v;
	vector<Target> r;
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
			false,
			ref_seqs.max_len(0, ref_seqs.size()),
			-1,
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
			r.emplace_back(block_id, ref_seqs[block_id], 0, nullptr);
		BlockId i = it.first->second;
		r[i].add_hit(hsp, hsp.begin());
	}

	return r;
}

}