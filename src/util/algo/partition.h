/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

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