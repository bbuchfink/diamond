/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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

template<typename _t>
struct FlatArray {

	typedef typename std::vector<_t>::iterator Iterator;
	typedef typename std::vector<_t>::const_iterator ConstIterator;

	FlatArray() {
		limits_.push_back(0);
	}

	void push_back(const _t &x) {
		data_.push_back(x);
		++limits_.back();
	}

	void push_back(ConstIterator begin, ConstIterator end) {
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

	void reserve(size_t n) {
		data_.reserve(n);
	}

	size_t size() const {
		return limits_.size() - 1;
	}

	size_t data_size() const {
		return data_.size();
	}

	ConstIterator begin(size_t i) const {
		return data_.cbegin() + limits_[i];
	}

	ConstIterator end(size_t i) const {
		return data_.cbegin() + limits_[i + 1];
	}

	ConstIterator cbegin(size_t i) const {
		return data_.cbegin() + limits_[i];
	}

	ConstIterator cend(size_t i) const {
		return data_.cbegin() + limits_[i + 1];
	}

	Iterator begin(size_t i) {
		return data_.begin() + limits_[i];
	}

	Iterator end(size_t i) {
		return data_.begin() + limits_[i + 1];
	}

	size_t count(size_t i) const {
		return limits_[i + 1] - limits_[i];
	}

private:

	std::vector<_t> data_;
	std::vector<size_t> limits_;

};