// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "basic/sequence.h"
#include <vector>

namespace DP {
struct FloatingScoring {
	const int* matrix; // Row-major, 32 letters per row, as in ScoreMatrix::matrix32().
	int gap_open, gap_extend;
};
// Score independent targets in SIMD lanes; no per-cell traceback storage.
std::vector<double> floating_scores(Sequence query, const std::vector<float>& bias,
	const std::vector<Sequence>& targets, FloatingScoring scoring);
}
