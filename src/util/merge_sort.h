/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef MERGE_SORT_H_
#define MERGE_SORT_H_

#include <algorithm>
#include <stddef.h>
#include "thread.h"

template<typename _it>
void merge_sort(_it begin, _it end, unsigned n_threads, unsigned level = 0)
{
	ptrdiff_t diff = end - begin;
	if(diff <= 1)
		return;

	if(1u << level >= n_threads) {
		std::sort(begin, end);
		return;
	}

	_it mid = begin + diff/2;
	thread *left = launch_thread(merge_sort<_it>, begin, mid, n_threads, level+1);
	thread *right = launch_thread(merge_sort<_it>, mid, end, n_threads, level+1);
	left->join();
	right->join();
	delete left;
	delete right;
	std::inplace_merge(begin, mid, end);
}

#endif /* MERGE_SORT_H_ */
