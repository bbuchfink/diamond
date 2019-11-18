/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

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

#include <vector>
#include <map>
#include "align.h"
#include "query_mapper.h"
#include "../data/reference.h"

using std::vector;
using std::map;

namespace ExtensionPipeline { namespace Swipe {

void Pipeline::run(Statistics &stat, const sequence *subjects, size_t subject_count)
{
	list<Hsp> hsp;
	if (subjects == nullptr) {
		vector<sequence> seqs;
		const size_t n = targets.size();
		seqs.reserve(n);
		for (size_t i = 0; i < n; ++i)
			seqs.push_back(ref_seqs::get()[targets[i].subject_block_id]);
		 hsp = DP::Swipe::swipe(query_seq(0), seqs.data(), seqs.data() + seqs.size(), raw_score_cutoff());
	}
	else
		hsp = DP::Swipe::swipe(query_seq(0), subjects, subjects + subject_count, raw_score_cutoff());	

	map<unsigned, unsigned> subject_idx;
	while (!hsp.empty()) {
		unsigned i;
		if (subjects == nullptr)
			i = hsp.begin()->swipe_target;
		else {
			const auto it = subject_idx.emplace(hsp.begin()->swipe_target, (unsigned)targets.size());
			if (it.second)
				targets.push_back(new Target(0, hsp.begin()->swipe_target));
			i = it.first->second;
		}
		targets[i].filter_score = hsp.begin()->score;
		list<Hsp> &l = targets[i].hsps;		
		l.splice(l.end(), hsp, hsp.begin());
	}
}

}}