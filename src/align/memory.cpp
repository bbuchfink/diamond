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

using std::pair;

namespace Extension {

Memory* memory;

Memory::Memory(size_t query_count):
	scores_(query_count * N, 0),
	count_(query_count, 0)
{
}

int& Memory::low_score(size_t query_id) {
	return scores_[query_id*N + (N - 1)];
}

int& Memory::mid_score(size_t query_id) {
	return scores_[query_id*N];
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
	const size_t cutoff = config.max_alignments, mid_size = (cutoff+1) / 2, low_size = cutoff - mid_size;
	size_t total_count = (size_t)count_[query_id], mid_count = std::min(mid_size, total_count), low_count = total_count - mid_count;
	int low_score = this->low_score(query_id), mid_score = this->mid_score(query_id);

	pair<size_t, int> r = update_range(begin, end, mid_size, mid_count, mid_score);
	this->mid_score(query_id) = mid_score;
	low_count = std::min(low_count + r.first, low_size);
	if (r.first >= low_size)
		low_score = r.second;
	
	update_range(begin, end, low_size, low_count, low_score);
	this->low_score(query_id) = low_score;
	this->count_[query_id] = mid_count + low_count;
}

}