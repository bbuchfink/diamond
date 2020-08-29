/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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
#include <algorithm>
#include "../basic/config.h"
#include "target.h"

using std::vector;

namespace Extension {

void rank_targets(vector<WorkTarget> &targets, double ratio, double factor)
{
	if (config.taxon_k && config.toppercent == 100.0)
		return;
	std::sort(targets.begin(), targets.end());
	if (targets.empty() || targets[0].filter_score == 0) {
		targets.clear();
		return;
	}

	int score = 0;
	if (config.toppercent < 100) {
		score = int((double)targets[0].filter_score * (1.0 - config.toppercent / 100.0) * ratio);
	}
	else {
		size_t min_idx = std::min(targets.size(), config.max_alignments);
		score = int((double)targets[min_idx - 1].filter_score * ratio);
	}
	score = std::max(score, 1);

	const size_t cap = (config.toppercent < 100 || config.max_alignments == std::numeric_limits<size_t>::max()) ? std::numeric_limits<size_t>::max() : size_t(config.max_alignments*factor);
	size_t i = 0;
	for (; i < targets.size(); ++i)
		if (targets[i].filter_score < score || i >= cap)
			break;

	if (config.benchmark_ranking)
		for (size_t j = i; j < targets.size(); ++j)
			targets[j].outranked = true;
	else
		targets.erase(targets.begin() + i, targets.end());
}

}