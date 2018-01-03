/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef TARGET_ITERATOR_H_
#define TARGET_ITERATOR_H_

#include "../dp.h"

template<int _n>
struct TargetIterator
{
	TargetIterator(vector<DpTarget>::const_iterator subject_begin, vector<DpTarget>::const_iterator subject_end) :
		next(0),
		n_targets(int(subject_end - subject_begin)),
		subject_begin(subject_begin)
	{
		for (; next < std::min(_n, n_targets); ++next) {
			pos[next] = 0;
			target[next] = next;
			active.push_back(next);
		}
	}
	char operator[](int i) const
	{
		return subject_begin[target[i]].seq[pos[i]];
	}
	__m128i get() const
	{
		char s[16];
		for (int i = 0; i < active.size(); ++i) {
			const int j = active[i];
			s[j] = (*this)[j];
		}
		return _mm_loadu_si128((const __m128i*)s);
	}
	bool init_target(int i, int j)
	{
		if (next < n_targets) {
			pos[j] = 0;
			target[j] = next++;
			return true;
		}
		active.erase(i);
		return false;
	}
	bool inc(int i)
	{
		++pos[i];
		if (pos[i] >= (int)subject_begin[target[i]].seq.length())
			return false;
		return true;
	}
	int pos[_n], target[_n], next, n_targets;
	Static_vector<int, _n> active;
	const vector<DpTarget>::const_iterator subject_begin;
};

#endif