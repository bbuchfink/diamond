/****
DIAMOND protein aligner
Copyright (C) 2020-2022 Max Planck Society for the Advancement of Science e.V.

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
#include <vector>
#include <stdint.h>
#define _REENTRANT
#include "../../lib/ips4o/ips4o.hpp"
#include "../algo/sort_helper.h"
#include "../algo/transform_iterator.h"
#include "../util.h"

template<typename T, typename I = int64_t>
struct FlatArray {

	typedef typename std::vector<T>::iterator DataIterator;
	typedef typename std::vector<T>::const_iterator DataConstIterator;

	FlatArray() {
		limits_.push_back(0);
	}

	FlatArray(std::vector<I>&& sizes):
		data_(sum(sizes)),
		limits_(sizes)
	{
	}

	FlatArray(std::vector<I>&& limits, std::vector<T>&& data) :
		data_(std::move(data)),
		limits_(std::move(limits))
	{}

	void push_back(const T &x) {
		data_.push_back(x);
		++limits_.back();
	}

	template<typename It>
	void push_back(const It begin, const It end) {
		data_.insert(data_.end(), begin, end);
		limits_.push_back(limits_.back() + (end - begin));
	}

	void next() {
		limits_.push_back(limits_.back());
	}

	void pop_back() {
		limits_.pop_back();
	}

	void clear() {
		data_.clear();
		limits_.clear();
		limits_.push_back(0);
	}

	void shrink_to_fit() {
		data_.shrink_to_fit();
		limits_.shrink_to_fit();
	}

	I size() const {
		return limits_.size() - 1;
	}

	I data_size() const {
		return data_.size();
	}

	DataConstIterator begin(I i) const {
		return data_.cbegin() + limits_[i];
	}

	DataConstIterator end(I i) const {
		return data_.cbegin() + limits_[i + 1];
	}

	DataConstIterator cbegin(I i) const {
		return data_.cbegin() + limits_[i];
	}

	DataConstIterator cend(I i) const {
		return data_.cbegin() + limits_[i + 1];
	}

	DataIterator begin(I i) {
		return data_.begin() + limits_[i];
	}

	DataIterator end(I i) {
		return data_.begin() + limits_[i + 1];
	}

	I count(I i) const {
		return limits_[i + 1] - limits_[i];
	}

	struct ConstIterator {
		ConstIterator(typename std::vector<size_t>::const_iterator limits, typename std::vector<T>::const_iterator data_begin):
			limits_(limits),
			data_begin_(data_begin)
		{}
		DataConstIterator begin(size_t i) const {
			return data_begin_ + limits_[i];
		}
		DataConstIterator end(size_t i) const {
			return data_begin_ + limits_[i + 1];
		}
		int64_t operator-(const ConstIterator& i) const {
			return limits_ - i.limits_;
		}
	private:
		typename std::vector<I>::const_iterator limits_;
		typename std::vector<T>::const_iterator data_begin_;
	};

	ConstIterator cbegin() const {
		return ConstIterator(limits_.cbegin(), data_.cbegin());
	}

	ConstIterator cend() const {
		return ConstIterator(limits_.cend() - 1, data_.cbegin());
	}

	struct Iterator {
		Iterator() {}
		Iterator(typename std::vector<I>::const_iterator limits, typename std::vector<T>::iterator data_begin) :
			limits_(limits),
			data_begin_(data_begin)
		{}
		DataIterator begin(I i) const {
			return data_begin_ + limits_[i];
		}
		DataIterator end(I i) const {
			return data_begin_ + limits_[i + 1];
		}
		int64_t operator-(const Iterator& i) const {
			return limits_ - i.limits_;
		}
	private:
		typename std::vector<I>::const_iterator limits_;
		typename std::vector<T>::iterator data_begin_;
	};

	Iterator begin() {
		return Iterator(limits_.cbegin(), data_.begin());
	}

	Iterator end() {
		return Iterator(limits_.cend() - 1, data_.begin());
	}

	void reserve(const I size, const I data_size) {
		data_.reserve(data_size);
		limits_.reserve(size + 1);
	}

	struct Range {
		Range(DataConstIterator begin, DataConstIterator end) :
			begin_(begin),
			end_(end)
		{}
		DataConstIterator begin() const {
			return begin_;
		}
		DataConstIterator end() const {
			return end_;
		}
	private:
		DataConstIterator begin_, end_;
	};

	Range operator[](I i) const {
		return Range { cbegin(i), cend(i) };
	}

	I max_count() const {
		I m = 0;
		for (I i = 0; i < size(); ++i)
			m = std::max(m, count(i));
		return m;
	}

private:

	I sum(std::vector<I>& sizes) {
		assert(sizes.back() == 0);
		I p = 0;
		for (typename std::vector<I>::iterator i = sizes.begin(); i < sizes.end(); ++i) {
			const I size = *i;
			*i = p;
			p += size;
		}
		return p;
	}

	std::vector<T> data_;
	std::vector<I> limits_;

};

template<typename It>
std::pair<FlatArray<typename It::value_type::second_type>, std::vector<typename It::value_type::first_type>> make_flat_array(const It begin, const It end, int num_threads)
{
	using Key = typename It::value_type::first_type;
	using Value = typename It::value_type::second_type;
#if _MSC_FULL_VER == 191627045 || !defined(NDEBUG)
	std::sort(begin, end);
#else
	ips4o::parallel::sort(begin, end, std::less<std::pair<Key, Value>>(), num_threads);
#endif
	auto it = merge_keys(begin, end, First<Key, Value>());
	int64_t n = 0;
	while (it.good()) {
		++n;
		++it;
	}
	std::pair<FlatArray<Value>, std::vector<Key>> r;
	r.first.reserve(n, end - begin);
	r.second.reserve(n);
	auto it2 = merge_keys(begin, end, First<Key, Value>());
	while (it2.good()) {
		r.second.push_back(it2.key());
		r.first.push_back(transform(it2.begin(), Second<Key, Value>()), transform(it2.end(), Second<Key, Value>()));
		++it2;
	}
	return r;
}


template<typename It>
FlatArray<typename It::value_type::second_type> make_flat_array_dense(const It begin, const It end, int num_threads)
{
	using Key = typename It::value_type::first_type;
	using Value = typename It::value_type::second_type;
#if _MSC_FULL_VER == 191627045 || !defined(NDEBUG)
	std:sort(begin, end);
#else
	ips4o::parallel::sort(begin, end, std::less<std::pair<Key, Value>>(), num_threads);
#endif
	const Key max_key = (end - 1)->first;
	FlatArray<Value> r;
	r.reserve(max_key + 1, end - begin);
	auto it = merge_keys(begin, end, First<Key, Value>());
	Key k = 0;
	while (it.good()) {
		while (k != it.key()) {
			r.next();
			++k;
		}
		r.push_back(transform(it.begin(), Second<Key, Value>()), transform(it.end(), Second<Key, Value>()));
		++k;
		++it;
	}
	return r;
}

template<typename T, typename GetKey>
FlatArray<T> make_flat_array_dense(std::vector<T>&& data, const typename T::Key key_end, int num_threads, GetKey get_key)
{
	std::vector<int64_t> limits;
	limits.push_back(0);
#if _MSC_FULL_VER == 191627045 || !defined(NDEBUG)
	std::sort(data.begin(), data.end());
#else
	ips4o::parallel::sort(data.begin(), data.end(), std::less<T>(), num_threads);
#endif
	limits.reserve(key_end + 1);
	auto it = merge_keys(data.cbegin(), data.cend(), get_key);
	typename T::Key k = 0;
	while (it.good()) {
		while (k != it.key()) {
			limits.push_back(limits.back());
			++k;
		}
		limits.push_back(limits.back() + it.count());
		++k;
		++it;
	}
	while (k++ < key_end)
		limits.push_back(limits.back());
	return FlatArray<T>(std::move(limits), std::move(data));
}