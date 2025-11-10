/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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