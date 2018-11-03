/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef COMPACT_ARRAY_H_
#define COMPACT_ARRAY_H_

#include <vector>
#include <limits.h>
#include "../io/deserializer.h"

using std::vector;

template <typename _t>
struct CompactArray
{

	CompactArray(Deserializer &in, size_t size, size_t data_size) :
		data_(data_size)
	{
		in.read(data_.data(), data_size);
		limits_.reserve(size + 1);
		limits_.push_back(0);
		Deserializer d(data_.data(), data_.data() + data_.size(), Deserializer::VARINT);
		_t x;
		for (size_t i = 0; i < size; ++i) {
			d >> x;
			size_t offset = d.data() - data_.data();
			if (offset > (size_t)UINT_MAX)
				throw std::runtime_error("Array size overflow.");
			limits_.push_back((unsigned)offset);
		}
		if (limits_.back() != data_size)
			throw std::runtime_error("Error loading CompactArray.");
	}

	_t operator[](size_t i) const
	{
		_t r;
		Deserializer(&data_[limits_[i]], &data_[limits_[i + 1]], Deserializer::VARINT) >> r;
		return r;
	}

	size_t size() const
	{
		return limits_.size() - 1;
	}

private:

	vector<char> data_;
	vector<unsigned> limits_;

};

#endif