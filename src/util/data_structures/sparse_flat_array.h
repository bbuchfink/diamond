/****
DIAMOND protein aligner
Copyright (C) 2022 Max Planck Society for the Advancement of Science e.V.

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
#include <unordered_map>
#include <tuple>
#include "flat_array.h"
#define _REENTRANT
#include "../../lib/ips4o/ips4o.hpp"
#include "../util.h"
#include "../algo/sort_helper.h"
#include "../algo/transform_iterator.h"

template<typename Key, typename T>
struct SparseFlatArray {

	using DataConstIterator = typename FlatArray<T>::DataConstIterator;
	using MapIt = typename std::unordered_map<Key, int64_t>::const_iterator;

	template<typename It>
	SparseFlatArray(const It begin, const It end, int num_threads)
	{
		ips4o::parallel::sort(begin, end, std::less<std::pair<Key, T>>(), num_threads);
		auto it = merge_keys(begin, end, First<Key, T>());
		int64_t n = 0;
		while (it.good()) {
			++n;
			++it;
		}
		map_.reserve(n);
		data_.reserve(n, end - begin);
		auto it2 = merge_keys(begin, end, First<Key, T>());
		while (it2.good()) {
			map_[it2.key()] = data_.size();
			data_.push_back(transform(it2.begin(), Second<Key, T>()), transform(it2.end(), Second<Key, T>()));
			++it2;
		}
	}

	bool empty() const {
		return map_.empty();
	}

	int64_t size() const {
		return map_.size();
	}

	int64_t data_size() const {
		return data_.data_size();
	}

	struct Iterator {
		Iterator(typename std::unordered_map<Key, int64_t>::iterator it, FlatArray<T>& data) :
			it_(it),
			data_(&data)
		{}
		Iterator operator++(int) {
			Iterator it = *this;
			++it_;
			return it;
		}
		bool operator==(const Iterator& it) const {
			return it_ == it.it_;
		}
		bool operator!=(const Iterator& it) const {
			return it_ != it.it_;
		}
		std::tuple<Key, typename FlatArray<T>::DataIterator, typename FlatArray<T>::DataIterator> operator*() const {
			return { it_->first, data_->begin(it_->second), data_->end(it_->second) };
		}
	private:
		typename std::unordered_map<Key, int64_t>::iterator it_;
		FlatArray<T>* data_;
	};

	struct ConstIterator {
		ConstIterator(MapIt it, const FlatArray<T>& data):
			it_(it),
			data_(&data)
		{}
		ConstIterator operator++(int) {
			Iterator it = *this;
			++it_;
			return it;
		}
		bool operator==(const ConstIterator& it) const {
			return it_ == it.it_;
		}
		bool operator!=(const ConstIterator& it) const {
			return it_ != it.it_;
		}
		std::tuple<Key, typename FlatArray<T>::DataConstIterator, typename FlatArray<T>::DataConstIterator> operator*() const {
			return { it_->first, data_->cbegin(it_->second), data_->cend(it_->second) };
		}
	private:
		typename std::unordered_map<Key, int64_t>::const_iterator it_;
		const FlatArray<T>* data_;
	};

	struct FlatIterator {
		FlatIterator(MapIt map_it, MapIt map_end, const FlatArray<T>& data) :
			map_it_(map_it),
			map_end_(map_end),
			data_it_(map_it != map_end ? data.cbegin(map_it->second) : DataConstIterator()),
			data_(&data)
		{}
		std::pair<Key, T> operator*() const {
			return std::make_pair(map_it_->first, *data_it_);
		}
		FlatIterator& operator++() {
			++data_it_;
			if (data_it_ >= data_->cend(map_it_->second)) {
				++map_it_;
				if(map_it_ != map_end_)
					data_it_ = data_->cbegin(map_it_->second);
			}
			return *this;
		}
		FlatIterator operator++(int) {
			FlatIterator it = *this;
			operator++();
			return it;
		}
		bool operator==(const FlatIterator& it) const {
			return map_it_ == it.map_it_ && (data_it_ == it.data_it_ || map_it_ == map_end_);
		}
		bool operator!=(const FlatIterator& it) const {
			return !(it == *this);
		}
	private:
		typename std::unordered_map<Key, int64_t>::const_iterator map_it_, map_end_;
		DataConstIterator data_it_;
		const FlatArray<T>* data_;
	};
	
	Iterator begin() {
		return Iterator(map_.begin(), data_);
	}

	Iterator end() {
		return Iterator(map_.end(), data_);
	}

	ConstIterator cbegin() const {
		return ConstIterator(map_.cbegin(), data_);
	}

	ConstIterator cend() const {
		return ConstIterator(map_.cend(), data_);
	}

	DataConstIterator cbegin(Key key) const {
		return data_.cbegin(map_.at(key));
	}

	DataConstIterator cend(Key key) const {
		return data_.cend(map_.at(key));
	}

	FlatIterator flat_begin() const {
		return FlatIterator(map_.cbegin(), map_.cend(), data_);
	}

	FlatIterator flat_end() const {
		return FlatIterator(map_.cend(), map_.cend(), data_);
	}

private:

	std::unordered_map<Key, int64_t> map_;
	FlatArray<T> data_;

};