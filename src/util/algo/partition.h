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
#include <stddef.h>

template<typename T = size_t>
struct Partition
{
	Partition() : parts(0), items(0), size_(0), remainder(0)
	{ }
	Partition(T items, T parts) : parts(std::min(parts, items)), items(items)
	{
		if (this->parts > 0) {
			size_ = items / this->parts;
			remainder = items % this->parts;
		}
		else {
			size_ = 0;
			remainder = 0;
		}
	}
	T begin(T i) const
	{
		T b = std::min(i, remainder); return b * (size_ + 1) + (i - b) * size_;
	}
	T end(T i) const
	{
		return begin(i) + size(i);
	}
	T size(T i) const
	{
		return i < remainder ? (size_ + 1) : size_;
	}
	const T parts;
private:
	T items, size_, remainder;
};