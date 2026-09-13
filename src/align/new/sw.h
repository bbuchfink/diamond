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

#include <limits>
#include <stdexcept>
#include <vector>
#include "stats/score_matrix.h"

namespace ExtensionPipeline {

enum Trace : uint8_t {
	STOP,
	DIAGONAL,
	INSERTION,
	DELETION,
	GAP_OPEN,
	GAP_EXTEND
};

struct TraceOperation {
	EditOperation op;
	Letter subject_letter;
};

static size_t cell(size_t query_pos, size_t target_pos, size_t columns) {
	return query_pos * columns + target_pos;
}

// Full-matrix Smith-Waterman with affine gaps. This deliberately favors
// readability over memory use and performance; the new pipeline currently
// aligns one protein query/target pair at a time.
static bool smith_waterman(const Sequence& query, const Sequence& target, Hsp& hsp) {
	if (query.empty() || target.empty())
		return false;

	const size_t rows = static_cast<size_t>(query.length()) + 1;
	const size_t columns = static_cast<size_t>(target.length()) + 1;
	const size_t matrix_size = rows * columns;
	const int negative_infinity = std::numeric_limits<int>::min() / 4;
	const int gap_open = score_matrix.gap_open() + score_matrix.gap_extend();
	const int gap_extend = score_matrix.gap_extend();

	std::vector<int> score(matrix_size, 0);
	std::vector<int> insertion_score(matrix_size, negative_infinity);
	std::vector<int> deletion_score(matrix_size, negative_infinity);
	std::vector<uint8_t> score_trace(matrix_size, STOP);
	std::vector<uint8_t> insertion_trace(matrix_size, GAP_OPEN);
	std::vector<uint8_t> deletion_trace(matrix_size, GAP_OPEN);

	int best_score = 0;
	size_t best_query_pos = 0, best_target_pos = 0;
	for (size_t i = 1; i < rows; ++i) {
		for (size_t j = 1; j < columns; ++j) {
			const size_t current = cell(i, j, columns);
			const int insertion_open = score[cell(i - 1, j, columns)] - gap_open;
			const int insertion_extend = insertion_score[cell(i - 1, j, columns)] - gap_extend;
			if (insertion_open >= insertion_extend) {
				insertion_score[current] = insertion_open;
				insertion_trace[current] = GAP_OPEN;
			}
			else {
				insertion_score[current] = insertion_extend;
				insertion_trace[current] = GAP_EXTEND;
			}

			const int deletion_open = score[cell(i, j - 1, columns)] - gap_open;
			const int deletion_extend = deletion_score[cell(i, j - 1, columns)] - gap_extend;
			if (deletion_open >= deletion_extend) {
				deletion_score[current] = deletion_open;
				deletion_trace[current] = GAP_OPEN;
			}
			else {
				deletion_score[current] = deletion_extend;
				deletion_trace[current] = GAP_EXTEND;
			}

			const int diagonal = score[cell(i - 1, j - 1, columns)] + score_matrix(query[i - 1], target[j - 1]);
			int current_score = 0;
			uint8_t trace = STOP;
			if (diagonal > current_score) {
				current_score = diagonal;
				trace = DIAGONAL;
			}
			if (insertion_score[current] > current_score) {
				current_score = insertion_score[current];
				trace = INSERTION;
			}
			if (deletion_score[current] > current_score) {
				current_score = deletion_score[current];
				trace = DELETION;
			}
			score[current] = current_score;
			score_trace[current] = trace;

			if (current_score > best_score) {
				best_score = current_score;
				best_query_pos = i;
				best_target_pos = j;
			}
		}
	}

	if (best_score == 0)
		return false;

	hsp.clear();
	hsp.backtraced = true;
	hsp.score = best_score;
	hsp.frame = 0;
	hsp.query_range.end_ = static_cast<Loc>(best_query_pos);
	hsp.subject_range.end_ = static_cast<Loc>(best_target_pos);

	std::vector<TraceOperation> transcript;
	size_t i = best_query_pos, j = best_target_pos;
	uint8_t state = score_trace[cell(i, j, columns)];
	while (state != STOP) {
		const size_t current = cell(i, j, columns);
		switch (state) {
		case DIAGONAL:
			transcript.push_back({ query[i - 1] == target[j - 1] ? op_match : op_substitution, target[j - 1] });
			--i;
			--j;
			state = score_trace[cell(i, j, columns)];
			break;
		case INSERTION:
			transcript.push_back({ op_insertion, Letter() });
			state = insertion_trace[current] == GAP_OPEN ? score_trace[cell(i - 1, j, columns)] : INSERTION;
			--i;
			break;
		case DELETION:
			transcript.push_back({ op_deletion, target[j - 1] });
			state = deletion_trace[current] == GAP_OPEN ? score_trace[cell(i, j - 1, columns)] : DELETION;
			--j;
			break;
		default:
			throw std::runtime_error("Invalid Smith-Waterman traceback state.");
		}
	}

	hsp.query_range.begin_ = static_cast<Loc>(i);
	hsp.subject_range.begin_ = static_cast<Loc>(j);
	hsp.query_source_range = hsp.query_range;
	hsp.subject_source_range = hsp.subject_range;

	EditOperation previous = op_match;
	bool have_previous = false;
	size_t query_pos = static_cast<size_t>(hsp.query_range.begin_);
	for (auto it = transcript.rbegin(); it != transcript.rend(); ++it) {
		if (it->op == op_match || it->op == op_substitution) {
			const Letter query_letter = query[query_pos];
			if (it->op == op_match) {
				hsp.transcript.push_back(op_match);
				++hsp.identities;
				++hsp.positives;
			}
			else {
				hsp.transcript.push_back(op_substitution, it->subject_letter);
				++hsp.mismatches;
				if (score_matrix(query_letter, it->subject_letter) > 0)
					++hsp.positives;
			}
			++query_pos;
		}
		else {
			if (it->op == op_insertion)
				hsp.transcript.push_back(op_insertion);
			else
				hsp.transcript.push_back(op_deletion, it->subject_letter);
			++hsp.gaps;
			if (!have_previous || previous != it->op)
				++hsp.gap_openings;
			if (it->op == op_insertion)
				++query_pos;
		}
		++hsp.length;
		previous = it->op;
		have_previous = true;
	}
	hsp.transcript.push_terminator();
	hsp.evalue = score_matrix.evalue(hsp.score, query.length(), target.length());
	hsp.bit_score = score_matrix.bitscore(hsp.score);
	hsp.corrected_bit_score = score_matrix.bitscore_corrected(hsp.score, query.length(), target.length());
	return true;
}

}