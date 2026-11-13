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
#include <functional>
#include "basic/value.h"
#include "basic/statistics.h"
#include "extension.h"

namespace MultipassPipeline {

// Receives a computed extension along with its score (0 if none was found) and its extent in the query
// and the target.
using ExtensionCallback = std::function<void(const Extension& extension, Score score, Loc query_len, Loc target_len)>;

// Computes the extensions [begin, end) in one call of the anchored swipe kernel and passes each of them to
// the callback once it is done. Every extension is an independent local alignment of the query and the target
// of its pair, starting at the anchor and running into its own direction, and only its score and extent are
// computed. band_max has to be at least the band of every extension (Extension::band).
// The extensions are loaded into the kernel in the given order. As the kernel computes all of its channels over
// the widest band of the extensions loaded so far, they should be ordered by band and have similar bands.
// common_scores is the 32x32 score table that the compact score tables of the pairs are expanded into (see
// CompactScoreTable).
void anchored_swipe(const Extension* const* begin, const Extension* const* end, Loc band_max, const int8_t* common_scores, Statistics& stats, const ExtensionCallback& callback);

}
