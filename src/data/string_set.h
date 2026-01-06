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
#include <cstddef>
#include <assert.h>
#include <vector>
#include <algorithm>

template<typename T, char padding_char, size_t padding_len = 1lu>
struct StringSetBase
{

	using Length = Loc;
	using Id = BlockId;
	using Pos = int64_t;

	enum { PERIMETER_PADDING = 256 };
	static const char DELIMITER = padding_char;

	StringSetBase():
		data_ (PERIMETER_PADDING, padding_char)
	{
		limits_.push_back(PERIMETER_PADDING);
	}

	void finish_reserve()
	{
		data_.resize(raw_len() + PERIMETER_PADDING);
		std::fill(data_.begin() + raw_len(), data_.end(), padding_char);
	}

	void reserve(size_t n)
	{
		limits_.push_back(raw_len() + n + padding_len);
	}

	void reserve(size_t entries, size_t length) {
		limits_.reserve(entries + 1);
		data_.reserve(length + 2 * PERIMETER_PADDING + entries * padding_len);
	}

	void clear() {
		limits_.resize(1);
		data_.resize(PERIMETER_PADDING);
	}

	void shrink_to_fit() {
		limits_.shrink_to_fit();
		data_.shrink_to_fit();
	}

	template<typename It>
	void push_back(It begin, It end)
	{
		assert(begin <= end);
		limits_.push_back(raw_len() + (end - begin) + padding_len);
		data_.insert(data_.end(), begin, end);
		data_.insert(data_.end(), padding_len, padding_char);
	}

	void append(const StringSetBase& s) {
		const Id n = s.size();
		if (n == 0)
			return;
		auto it = s.limits_.cbegin() + 1;
		assert(raw_len() >= s.limits_.front());
		const Pos offset = raw_len() - s.limits_.front();
		//limits_.reserve(limits_.size() + n);
		for (Id i = 0; i < n; ++i)
			limits_.push_back(*it++ + offset);
		data_.insert(data_.end(), s.ptr(0), s.end(n - 1) + 1);
	}

	template<typename It>
	void assign(const size_t i, const It begin, const It end) {
		std::copy(begin, end, ptr(i));
		std::fill(ptr(i) + (end - begin), ptr(i) + (end - begin) + padding_len, padding_char);
	}

	void fill(size_t n, T v)
	{
		limits_.push_back(raw_len() + n + padding_len);
		data_.insert(data_.end(), n, v);
		data_.insert(data_.end(), padding_len, padding_char);
	}

	T* ptr(size_t i)
	{ return &data_[limits_[i]]; }

	const T* ptr(size_t i) const
	{ return &data_[limits_[i]]; }

	const T* end(size_t i) const {
		return &data_[limits_[i + 1] - padding_len];
	}

	size_t check_idx(size_t i) const
	{
		if (limits_.size() < i + 2)
			throw std::runtime_error("Sequence set index out of bounds.");
		return i;
	}

	Length length(size_t i) const
	{
		return Length(limits_[i + 1] - limits_[i] - padding_len);
	}

	Id size() const
	{ return Id(limits_.size() - 1); }

	bool empty() const {
		return limits_.size() <= 1;
	}

	int64_t raw_len() const
	{ return limits_.back(); }

	int64_t mem_size() const {
		return data_.size() * sizeof(T) + limits_.size() * sizeof(int64_t);
	}

	int64_t letters() const
	{ return raw_len() - size() - PERIMETER_PADDING; }

	T* data(uint64_t p = 0)
	{ return &data_[p]; }

	const T* data(uint64_t p = 0) const
	{ return &data_[p]; }

	size_t position(const T* p) const
	{ return p - data(); }

	Pos position(Id i, Length j) const
	{ return limits_[i] + j; }

	std::pair<Id, Length> local_position(int64_t p) const
	{
		auto i = std::upper_bound(limits_.begin(), limits_.end(), p) - limits_.begin() - 1;
		return std::pair<Id, Length>(Id(i), Length(p - limits_[i]));
	}

	template<typename It, typename Out, typename Cmp>
	void local_position_batch(It begin, It end, Out out, Cmp cmp) const {
		batch_binary_search(begin, end, limits_.begin(), limits_.end(), out, cmp);
	}

	const T* operator[](size_t i) const
	{
		return ptr(i);
	}

	const T* back() const {
		return ptr(limits_.size() - 2);
	}

	typename std::vector<int64_t>::const_iterator limits_begin() const {
		return limits_.begin();
	}

	typename std::vector<int64_t>::const_iterator limits_end() const {
		return limits_.end();
	}

	struct ConstIterator {

		ConstIterator(const T* data, const int64_t* limits):
			data_(data),
			limits_(limits)
		{}

		using iterator_category = std::random_access_iterator_tag;
		using difference_type = ptrdiff_t;
		using value_type = T;
		using pointer = T*;
		using reference = T&;

		ptrdiff_t operator-(const ConstIterator& it) const {
			return limits_ - it.limits_;
		}

		bool operator==(const ConstIterator& it) const {
			return limits_ == it.limits_;
		}

		bool operator!=(const ConstIterator& it) const {
			return limits_ != it.limits_;
		}

		ConstIterator operator+(ptrdiff_t d) const {
			return { data_ + *(limits_ + d) - *limits_, limits_ + d };
		}

		bool operator<(const ConstIterator& it) const {
			return limits_ < it.limits_;
		}

		ConstIterator& operator+=(ptrdiff_t d) {
			data_ += *(limits_ + d) - *limits_;
			limits_ += d;
			return *this;
		}

		ConstIterator& operator++() {
			this->operator+=(1);
			return *this;
		}

		std::pair<const T*, int64_t> operator[](const ptrdiff_t i) const {
			return { data_ + limits_[i] - limits_[0], limits_[i + 1] - limits_[i] - padding_len };
		}

		std::pair<const T*, int64_t> operator*() const {
			return this->operator[](0);
		}

	private:

		const T* data_;
		const int64_t* limits_;

	};

	ConstIterator cbegin() const {
		return ConstIterator(ptr(0), limits_.data());
	}

	ConstIterator cend() const {
		return ConstIterator(nullptr, &limits_[size()]);
	}

	template<typename It>
	StringSetBase subset(It begin, It end) const {
		StringSetBase r;
		r.limits_.reserve(end - begin);
		for (It i = begin; i != end; ++i)
			r.reserve(length(*i));
		r.finish_reserve();
		Id n = 0;
		for (It i = begin; i != end; ++i, ++n) {
			r.assign(n, ptr(*i), this->end(*i));
		}
		return r;
	}

private:

	std::vector<T> data_;
	std::vector<Pos> limits_;

};

using StringSet = StringSetBase<char, '\0'>;