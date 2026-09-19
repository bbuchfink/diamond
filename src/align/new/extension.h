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

#pragma once
#include <algorithm>
#include <cfloat>
#include <functional>
#include <mutex>
#include <vector>
#include "basic/sequence.h"
#include "util/geo/interval.h"
#include "util/hsp/approx_hsp.h"

namespace ExtensionPipeline {

struct QueryState;

// Extension of an anchor to a local alignment of a query/target pair. The target task
// submitting it sets the input members, the anchored swipe sets the result members.
struct Extension {
	Extension() :
		state(nullptr),
		target_id(0),
		score(0),
		evalue(DBL_MAX)
	{}
	Extension(QueryState* state, BlockId target_id, Sequence query, Sequence target, const Anchor& anchor) :
		state(state),
		target_id(target_id),
		query(query),
		target(target),
		anchor(anchor),
		score(0),
		evalue(DBL_MAX)
	{}
	// Widest diagonal range of the anchor over both extension directions, including the anchor diagonal itself.
	// Unset ranges count as empty.
	Loc diag_spread() const {
		const Loc d = anchor.diag();
		return std::max(std::max(anchor.d_max_right, d) - std::min(anchor.d_min_right, d),
			std::max(anchor.d_max_left, d) - std::min(anchor.d_min_left, d));
	}
	QueryState* state;            // not used by the anchored swipe
	BlockId target_id;
	Sequence query, target;
	Anchor anchor;
	Score score;                  // 0 if no alignment passing the filters was found
	Interval query_range, target_range;
	double evalue;
};

// Queue of the extensions waiting to be computed. The anchored swipe pops all of them at
// once and holds the lock only while doing so, so the queue can be filled again while the
// extensions are being computed.
struct ExtensionQueue {
	std::mutex mtx;
	std::vector<Extension> extensions;
};

// Receives the extensions popped from the queue once they have been computed.
using ExtensionCallback = std::function<void(std::vector<Extension>&)>;

}