/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#include "align.h"
#include "query_mapper.h"
#include "../data/reference.h"

namespace ExtensionPipeline { namespace Swipe {

void Pipeline::run(Statistics &stat)
{
	const size_t n = targets.size();
	vector<DpTarget> seqs(n);
	for (size_t i = 0; i < n; ++i) {
		seqs[i].seq = ref_seqs::get()[targets[i].subject_id];
	}
	vector<int> scores(n);
	swipe(query_seq(0), seqs.begin(), seqs.end(), scores.begin());
	for (size_t i = 0; i < n; ++i)
		targets[i].hsps.push_back(Hsp(scores[i]));
}

}}