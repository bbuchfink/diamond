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

	using Length = Loc;
	using Id = BlockId;
	using Pos = int64_t;

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

	void append(const StringSetBase& s) {
		reserve(size() + s.size(), letters() + s.letters());
		for (Id i = 0; i < s.size(); ++i)
			push_back(s.ptr(i), s.end(i));
	}

	template<typename It>
	void assign(const size_t i, const It begin, const It end) {
		std::copy(begin, end, ptr(i));
		std::fill(ptr(i) + (end - begin), ptr(i) + (end - begin) + _padding, _pchar);
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

	Length length(size_t i) const
	{
		return Length(limits_[i + 1] - limits_[i] - _padding);
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
			return { data_ + limits_[i] - limits_[0], limits_[i + 1] - limits_[i] - _padding };
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