/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

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

template<typename T>
struct DoubleBuffer
{

	inline void init(size_t size, size_t padding, size_t padding_front, T init)
	{
		const size_t total = size + padding + padding_front;
		data_.clear();
		data_.resize(total * 2);
		ptr1 = &data_[padding_front];
		ptr2 = &data_[total + padding_front];
		for (size_t i = 0; i < total * 2; ++i)
			data_[i] = init;
	}

	inline std::pair<T*, T*> get(int)
	{
		std::swap(ptr1, ptr2); return std::pair<T*, T*>(ptr2, ptr1);
	}

	inline T* last()
	{
		return ptr1;
	}

private:
	T* ptr1, * ptr2;
	std::vector<T> data_;

};