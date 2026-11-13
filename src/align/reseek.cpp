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

#include <array>
#include <algorithm>
#include <functional>
#include <vector>
#include <stdexcept>
#include "basic/config.h"
#include "reseek.h"
#include "search/hamming/finger_print.h"
#include "util/simd/dispatch.h"

using std::vector;
using std::runtime_error;

namespace Extension { namespace DISPATCH_ARCH {

vector<SeedHit> reseek_diags(const Query& query, Sequence target, unsigned hamming_filter_id) {
	using ::DISPATCH_ARCH::FingerPrint;
	vector<SeedHit> hits;
	const Loc offset = target.length() - 1, t = query.sequence[0].length() - 1 + target.length();
	vector<Loc> last_hit(t, -1);
	vector<bool> double_hit(t, false);
	auto f = [&hits, &last_hit, &double_hit, offset, &query, target, hamming_filter_id](Loc query_loc, Loc target_loc) {
		const Loc d = query_loc - target_loc + offset;
		if (d < 0 || d >= (Loc)last_hit.size())
			throw runtime_error("Invalid diagonal");
		if (last_hit[d] != -1 && target_loc - last_hit[d] < query.seqindex->k()) {
			return;
		}
		if (last_hit[d] == -1 || target_loc - last_hit[d] > config.double_hit_window) {
			double_hit[d] = false;
		} else if (!double_hit[d]) {
			SeedHit hit;
			const Loc s = (target_loc - last_hit[d]) / 2;
			hit.i = std::max(query_loc - s, 0);
			hit.j = std::max(target_loc - s, 0);
			hit.frame = 0;
			hit.score = 0;
			if (hamming_filter_id == 0) {
				hits.push_back(hit);
			}
			else {
				std::array<char, 48> fq, fs;
				FingerPrint::load(query.sequence[0].data() + hit.i, &fq);
				FingerPrint::load(target.data() + hit.j, &fs);
				if (FingerPrint(fq).match(FingerPrint(fs)) >= hamming_filter_id)
					hits.push_back(hit);
			}
			double_hit[d] = true;
		}
		last_hit[d] = target_loc;
	};
	query.seqindex->scan(target, f);
	return hits;
}

}

DISPATCH_3(vector<SeedHit>, reseek_diags, const Query&, query, Sequence, target, unsigned, hamming_filter_id)

}
