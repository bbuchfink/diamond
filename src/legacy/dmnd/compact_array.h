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

#pragma once
#include <vector>
#include <limits.h>
#include "util/io/deserializer.h"
#include "util/algo/varint.h"

struct CompactArray
{

	CompactArray(Deserializer &in, size_t size, size_t data_size) :
		data_(data_size)
	{
		in.read(data_.data(), data_size);
		if(data_size > (size_t)UINT_MAX)
			init(size, data_size, limits64_);
		else
			init(size, data_size, limits_);
	}

	std::vector<int32_t> operator[](size_t i) const
	{
		return data_.size() > (size_t)UINT_MAX ? get(i, limits64_) : get(i, limits_);
	}

	size_t size() const
	{
		return (data_.size() > (size_t)UINT_MAX ? limits64_.size() : limits_.size()) - 1;
	}

private:

	template<typename Int>
	void init(int64_t size, int64_t data_size, std::vector<Int>& limits) {
		limits.reserve(size + 1);
		limits.push_back(0);
		const char* ptr = data_.data();
		for (int64_t i = 0; i < size; ++i) {
			ptr = skip_vec(ptr);
			const int64_t offset = ptr - data_.data();
			if (offset > (int64_t)std::numeric_limits<Int>::max())
				throw std::runtime_error("Array size overflow.");
			limits.push_back((Int)offset);
		}
		if ((int64_t)limits.back() != data_size)
			throw std::runtime_error("Error loading CompactArray.");
	}

	template<typename Int>
	std::vector<int32_t> get(size_t i, const std::vector<Int>& limits) const
	{
		return read_vec(&data_[limits[i]]);
	}

	static std::vector<int32_t> read_vec(const char* ptr) {
		uint32_t n;
		std::tie(n, ptr) = read_varuint32(ptr);
		std::vector<int32_t> out;
		out.reserve(n);
		for (uint32_t i = 0; i < n; ++i) {
			uint32_t v;
			std::tie(v, ptr) = read_varuint32(ptr);
			out.push_back(v);
		}
		return out;
	}

	static const char* skip_vec(const char* ptr) {
		uint32_t n;
		std::tie(n, ptr) = read_varuint32(ptr);
		for (uint32_t i = 0; i < n; ++i) {
			uint32_t v;
			std::tie(v, ptr) = read_varuint32(ptr);
		}
		return ptr;
	}

	std::vector<char> data_;
	std::vector<unsigned> limits_;
	std::vector<int64_t> limits64_;

};
