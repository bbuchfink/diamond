/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include <assert.h>
#include <vector>
#include <stddef.h>
#include "../basic/sequence.h"
#include "../util/algo/binary_search.h"

template<typename T, char _pchar, size_t _padding = 1lu>
struct StringSetBase
{

	enum { PERIMETER_PADDING = 256 };
	static const char DELIMITER = _pchar;

	StringSetBase():
		data_ (PERIMETER_PADDING, _pchar)
	{
		limits_.push_back(PERIMETER_PADDING);
	}

	void finish_reserve()
	{
		data_.resize(raw_len() + PERIMETER_PADDING);
		std::fill(data_.begin() + raw_len(), data_.end(), _pchar);
	}

	void reserve(size_t n)
	{
		limits_.push_back(raw_len() + n + _padding);
	}

	void reserve(size_t entries, size_t length) {
		limits_.reserve(entries + 1);
		data_.reserve(length + 2 * PERIMETER_PADDING + entries * _padding);
	}

	void clear() {
		limits_.resize(1);
		data_.resize(PERIMETER_PADDING);
	}

	void shrink_to_fit() {
		limits_.shrink_to_fit();
		data_.shrink_to_fit();
	}

	template<typename _it>
	void push_back(_it begin, _it end)
	{
		assert(begin <= end);
		limits_.push_back(raw_len() + (end - begin) + _padding);
		data_.insert(data_.end(), begin, end);
		data_.insert(data_.end(), _padding, _pchar);
	}

	template<typename It>
	void assign(const size_t i, const It begin, const It end) {
		std::copy(begin, end, ptr(i));
	}

	void fill(size_t n, T v)
	{
		limits_.push_back(raw_len() + n + _padding);
		data_.insert(data_.end(), n, v);
		data_.insert(data_.end(), _padding, _pchar);
	}

	T* ptr(size_t i)
	{ return &data_[limits_[i]]; }

	const T* ptr(size_t i) const
	{ return &data_[limits_[i]]; }

	const T* end(size_t i) const {
		return &data_[limits_[i + 1] - _padding];
	}

	size_t check_idx(size_t i) const
	{
		if (limits_.size() < i + 2)
			throw std::runtime_error("Sequence set index out of bounds.");
		return i;
	}

	size_t length(size_t i) const
	{ return limits_[i+1] - limits_[i] - _padding; }

	size_t size() const
	{ return limits_.size() - 1; }

	bool empty() const {
		return limits_.size() <= 1;
	}

	size_t raw_len() const
	{ return limits_.back(); }

	size_t letters() const
	{ return raw_len() - size() - PERIMETER_PADDING; }

	T* data(uint64_t p = 0)
	{ return &data_[p]; }

	const T* data(uint64_t p = 0) const
	{ return &data_[p]; }

	size_t position(const T* p) const
	{ return p - data(); }

	size_t position(size_t i, size_t j) const
	{ return limits_[i] + j; }

	std::pair<size_t, size_t> local_position(size_t p) const
	{
		size_t i = std::upper_bound(limits_.begin(), limits_.end(), p) - limits_.begin() - 1;
		return std::pair<size_t, size_t>(i, p - limits_[i]);
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

	typename std::vector<size_t>::const_iterator limits_begin() const {
		return limits_.begin();
	}

	typename std::vector<size_t>::const_iterator limits_end() const {
		return limits_.end();
	}

	struct ConstIterator {

		ConstIterator(const T* data, const size_t* limits):
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

		std::pair<const T*, size_t> operator[](const ptrdiff_t i) const {
			return { data_ + limits_[i] - limits_[0], limits_[i + 1] - limits_[i] - _padding };
		}

	private:

		const T* data_;
		const size_t* limits_;

	};

	ConstIterator cbegin() const {
		return ConstIterator(ptr(0), limits_.data());
	}

	ConstIterator cend() const {
		return ConstIterator(nullptr, &limits_[size()]);
	}

private:

	std::vector<T> data_;
	std::vector<size_t> limits_;

};

using StringSet = StringSetBase<char, '\0'>;