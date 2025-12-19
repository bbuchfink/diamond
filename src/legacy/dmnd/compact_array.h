/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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