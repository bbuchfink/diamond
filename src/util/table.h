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
#include <string>
#include <utility>
#include <ostream>
#include <iomanip>
#include <math.h>

namespace Util {

struct Table {

	Table():
		max_len_(0)
	{}

	Table& operator()(const std::string& s, const std::string& t) {
		data_.emplace_back(s, t);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, long long n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, unsigned long long n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, long n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, unsigned long n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, int n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, unsigned int n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, double n, const char* unit = "") {
		data_.emplace_back(s, (n >= 100.0 ? std::to_string((int64_t)round(n)) : std::to_string(n)) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	friend std::ostream& operator<<(std::ostream& str, const Table& table) {
		for (const auto& l : table.data_)
			str << std::setw(table.max_len_) << l.first << "  " << l.second << std::endl;
		return str;
	}

private:

	std::vector<std::pair<std::string, std::string>> data_;
	size_t max_len_;

};

}