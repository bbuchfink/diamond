// SPDX-License-Identifier: GPL-3.0-or-later
#include "dp/scalar/floating_sw.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <string>

static void require(bool condition, const char* message) {
	if (!condition) throw std::runtime_error(message);
}

int main() {
	try {
		const auto fractional = DP::floating_sw(2, 2, 11, 1,
			[](size_t, size_t) { return 2.0 - 0.49; });
		require(std::abs(fractional.score - 3.02) < 1e-12, "Fractional corrections were lost.");
		auto substitution = [](size_t i, size_t j) { return i == j ? 2.0 : (i == 0 && j == 3 ? 7.0 : -20.0); };
		const auto unrounded_path = DP::floating_sw(4, 4, 11, 1,
			[&](size_t i, size_t j) { return substitution(i, j) - 0.49; });
		const auto rounded_path = DP::floating_sw(4, 4, 11, 1, substitution);
		require(std::abs(unrounded_path.score - 6.51) < 1e-12
			&& unrounded_path.operations.size() == 1 && rounded_path.operations.size() == 4,
			"Fractional scores must select their own optimal alignment, not rescore the integer path.");
		for (bool reverse : { false, true }) {
			const std::string q = reverse ? "AWQQAA" : "AWAA", t = reverse ? "AWAA" : "AWQQAA";
			const auto gapped = DP::floating_sw(q.size(), t.size(), 1.0, 0.25,
				[&](size_t i, size_t j) { return q[i] == t[j] ? 3.25 : -9.0; });
			require(std::abs(gapped.score - 11.5) < 1e-12, "Affine gap cost is incorrect.");
			require(std::string(gapped.operations.begin(), gapped.operations.end()) == (reverse ? "MMIIMM" : "MMDDMM"),
				"Gap traceback is incorrect.");
		}
		const auto empty = DP::floating_sw(0, 4, 1, 1, [](size_t, size_t) { return 1.0; });
		require(empty.score == 0 && empty.operations.empty(), "Empty alignment is incorrect.");
		const auto negative = DP::floating_sw(2, 2, 1, 1, [](size_t, size_t) { return -0.1; });
		require(negative.score == 0 && negative.operations.empty(), "Local reset is incorrect.");
		// Exhaustively enumerate all positive local paths on small random matrices.
		// This is independent of the dynamic-programming recurrence and also checks
		// that the recorded traceback reproduces the selected fractional score.
		std::mt19937 rng(19);
		for (int trial = 0; trial < 200; ++trial) {
			const size_t n = 1 + rng() % 5, m = 1 + rng() % 5;
			std::vector<double> substitutions(n * m);
			for (auto& x : substitutions) x = (int(rng() % 37) - 18) / 7.0;
			const double open = 1.25, extend = 0.4;
			auto sub = [&](size_t i, size_t j) { return substitutions[i * m + j]; };
			double best = 0;
			std::function<void(size_t, size_t, double, char)> visit = [&](size_t i, size_t j, double s, char previous) {
				if (s <= 0) return;
				best = std::max(best, s);
				if (i < n && j < m) visit(i + 1, j + 1, s + sub(i, j), 'M');
				if (i < n) visit(i + 1, j, s - extend - (previous == 'I' ? 0 : open), 'I');
				if (j < m) visit(i, j + 1, s - extend - (previous == 'D' ? 0 : open), 'D');
			};
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < m; ++j) visit(i + 1, j + 1, sub(i, j), 'M');
			const auto a = DP::floating_sw(n, m, open, extend, sub);
			require(std::abs(a.score - best) < 1e-10, "Score differs from exhaustive path enumeration.");
			size_t i = a.query_begin, j = a.target_begin;
			double score = 0;
			char previous = 'M';
			for (char op : a.operations) {
				if (op == 'M') { score += sub(i, j); ++i; ++j; }
				else {
					score -= extend + (previous == op ? 0 : open);
					if (op == 'I') ++i; else ++j;
				}
				previous = op;
			}
			require(std::abs(score - a.score) < 1e-10 && i == a.query_end && j == a.target_end,
				"Traceback does not reproduce the score and coordinates.");
		}
	}
	catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
	return 0;
}
