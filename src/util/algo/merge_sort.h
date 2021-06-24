/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#pragma once
#include <algorithm>
#include <stddef.h>
#include <thread>
#include <iterator>

template<typename It, typename Cmp = std::less<typename std::iterator_traits<It>::value_type>>
void merge_sort(It begin, It end, unsigned n_threads, const Cmp& cmp = std::less<typename std::iterator_traits<It>::value_type>(), unsigned level = 0)
{
	ptrdiff_t diff = end - begin;
	if(diff <= 1)
		return;

	if(1u << level >= n_threads) {
		std::sort(begin, end, cmp);
		return;
	}

	It mid = begin + diff/2;
	std::thread *left = new std::thread(merge_sort<It, Cmp>, begin, mid, n_threads, cmp, level+1);
	std::thread *right = new std::thread(merge_sort<It, Cmp>, mid, end, n_threads, cmp, level+1);
	left->join();
	right->join();
	delete left;
	delete right;
	std::inplace_merge(begin, mid, end, cmp);
}
