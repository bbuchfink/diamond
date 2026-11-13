// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace DP {

struct FloatingAlignment {
	double score = 0.0;
	size_t query_begin = 0, query_end = 0, target_begin = 0, target_end = 0;
	std::vector<char> operations;
};

// Scalar reference implementation. Gap cost is open + length * extend,
// matching DIAMOND's integer kernels. Traceback stores decisions rather than
// comparing floating-point sums for equality.
template<typename Substitution>
FloatingAlignment floating_sw(size_t query_len, size_t target_len,
	double gap_open, double gap_extend, Substitution substitution) {
	const size_t width = target_len + 1;
	if (width == 0 || query_len == std::numeric_limits<size_t>::max()
		|| query_len + 1 > std::vector<uint8_t>().max_size() / width)
		throw std::length_error("Floating-point alignment matrix is too large.");
	std::vector<uint8_t> trace((query_len + 1) * width, 0);
	std::vector<double> previous(width, 0.0), current(width, 0.0);
	std::vector<double> vertical(width, -std::numeric_limits<double>::infinity());
	FloatingAlignment out;
	for (size_t i = 1; i <= query_len; ++i) {
		current[0] = 0.0;
		double horizontal = -std::numeric_limits<double>::infinity();
		for (size_t j = 1; j <= target_len; ++j) {
			uint8_t flags = 0;
			const double h_open = current[j - 1] - gap_open - gap_extend;
			const double h_extend = horizontal - gap_extend;
			if (h_extend > h_open) flags |= 4;
			horizontal = std::max(h_open, h_extend);
			const double v_open = previous[j] - gap_open - gap_extend;
			const double v_extend = vertical[j] - gap_extend;
			if (v_extend > v_open) flags |= 8;
			vertical[j] = std::max(v_open, v_extend);
			double score = 0.0;
			uint8_t direction = 0;
			const double diagonal = previous[j - 1] + substitution(i - 1, j - 1);
			if (diagonal > score) { score = diagonal; direction = 1; }
			if (horizontal > score) { score = horizontal; direction = 2; }
			if (vertical[j] > score) { score = vertical[j]; direction = 3; }
			current[j] = score;
			trace[i * width + j] = flags | direction;
			if (score > out.score) {
				out.score = score;
				out.query_end = i;
				out.target_end = j;
			}
		}
		previous.swap(current);
	}
	size_t i = out.query_end, j = out.target_end;
	uint8_t state = 0;
	while (i > 0 && j > 0) {
		const uint8_t flags = trace[i * width + j];
		if (state == 0) {
			const uint8_t direction = flags & 3;
			if (direction == 0) break;
			if (direction == 1) {
				out.operations.push_back('M');
				--i; --j;
				continue;
			}
			state = direction;
		}
		if (state == 2) {
			out.operations.push_back('D');
			--j;
			state = (flags & 4) ? 2 : 0;
		}
		else {
			out.operations.push_back('I');
			--i;
			state = (flags & 8) ? 3 : 0;
		}
	}
	out.query_begin = i;
	out.target_begin = j;
	std::reverse(out.operations.begin(), out.operations.end());
	return out;
}

}
