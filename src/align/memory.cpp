/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include <utility>
#include "../basic/config.h"
#include "target.h"
#include "../util/util.h"

using std::pair;

namespace Extension {

Memory* memory = nullptr;

Memory::Memory(size_t query_count):
	N(config.memory_intervals),
	scores_(query_count * N, 0),
	count_(query_count, 0),
	ranking_low_score_(query_count, 0),
	ranking_failed_count_(query_count, 0)
{
}

int& Memory::low_score(size_t query_id) {
	return scores_[query_id*N + (N - 1)];
}

int& Memory::mid_score(size_t query_id) {
	return scores_[query_id*N];
}

int& Memory::min_score(size_t query_id, size_t i) {
	return scores_[query_id*N + i];
}

size_t Memory::count(size_t query_id) const {
	return (size_t)count_[query_id];
}

static pair<size_t, int> update_range(vector<Target>::const_iterator& begin, vector<Target>::const_iterator end, size_t size, size_t& count, int& low_score) {
	if (begin >= end)
		return { 0, 0 };
	auto it = begin;
	size_t n = 0, cut = 0, total = count;
	const int low = low_score;
	while (it < end && (it->filter_score > low || total < size) && n < size) {
		++it;
		++n;
		if (total < size)
			++total;
		else
			++cut;
	}
	if (n == 0)
		return { 0, 0 };
	begin = it;
	--it;
	if (n == total)
		low_score = it->filter_score;
	else if (count < total)
		low_score = std::min(low_score, it->filter_score);
	count = total;
	return { cut, low };
}

void Memory::update(size_t query_id, std::vector<Target>::const_iterator begin, std::vector<Target>::const_iterator end) {
	if (config.no_query_memory)
		return;

	const size_t cutoff = config.max_alignments;
	partition<size_t> p(cutoff, N);
	size_t total = count_[query_id], overflow_count = 0;
	int overflow_score = 0;
	for (size_t i = 0; i < p.parts; ++i) {
		const size_t size = p.getCount(i);
		size_t count = std::min(total, size);
		int low_score = this->min_score(query_id, i);
		if (overflow_count >= size)
			low_score = std::max(low_score, overflow_score);
		auto it = begin;
		pair<size_t, int> r = update_range(begin, end, size, count, low_score);
		
		overflow_count = r.first;
		overflow_score = r.second;
		this->min_score(query_id, i) = low_score;
		const size_t n = begin - it;
		total += n;
		total -= count;
		count_[query_id] += n;
	}

	count_[query_id] = std::min(count_[query_id], (int)cutoff);
}

void Memory::update_failed_count(size_t query_id, size_t failed_count, int ranking_low_score) {
	if (ranking_low_score >= ranking_low_score_[query_id]) {
		ranking_low_score_[query_id] = ranking_low_score;
		ranking_failed_count_[query_id] = failed_count;
	}
}


}