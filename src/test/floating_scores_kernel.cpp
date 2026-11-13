// SPDX-License-Identifier: GPL-3.0-or-later
#include "dp/floating_scores.h"
#include "dp/scalar/floating_sw.h"
#include <array>
#include <cmath>
#include <iostream>
#include <random>

// The score kernel is also compiled in isolation, with a deterministic test
// matrix and without runtime dispatch or global application configuration.
namespace DP { namespace FLOAT_TEST {
std::vector<double> floating_scores(Sequence query, const std::vector<float>& bias,
	const std::vector<Sequence>& targets, FloatingScoring scoring);
}}

int main() {
	try {
		std::mt19937 rng(1003);
		std::array<int, 32 * 32> matrix;
		for (auto& score : matrix) score = int(rng() % 17) - 8;
		for (int trial = 0; trial < 200; ++trial) {
			std::vector<Letter> query(rng() % 32);
			std::vector<float> bias(query.size());
			for (size_t i = 0; i < query.size(); ++i) {
				query[i] = Letter(rng() % 26);
				bias[i] = float(int(rng() % 61) - 30) / 13.0f;
			}
			std::vector<std::vector<Letter>> data(rng() % 18);
			std::vector<Sequence> targets;
			for (auto& target : data) {
				target.resize(rng() % 40);
				for (auto& letter : target) letter = Letter(rng() % 26);
				targets.emplace_back(target);
			}
			const int open = rng() % 8, extend = 1 + rng() % 3;
			const auto scores = DP::FLOAT_TEST::floating_scores(Sequence(query), bias, targets,
				{ matrix.data(), open, extend });
			for (size_t k = 0; k < targets.size(); ++k) {
				const auto reference = DP::floating_sw(query.size(), targets[k].length(), open, extend,
					[&](size_t i, size_t j) { return double(matrix[size_t(query[i]) * 32 + targets[k][j]]) + double(bias[i]); });
				if (scores[k] != reference.score)
					throw std::runtime_error("SIMD score differs from scalar reference.");
			}
		}
	}
	catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
	return 0;
}
