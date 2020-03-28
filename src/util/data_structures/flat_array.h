/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

#ifndef FLAT_ARRAY_H_
#define FLAT_ARRAY_H_

#include <vector>

template<typename _t>
struct FlatArray {

	FlatArray() {
		limits_.push_back(0);
	}

	void push_back(const _t &x) {
		data_.push_back(x);
		++limits_.back();
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

	size_t size() const {
		return limits_.size() - 1;
	}

	const _t* begin(size_t i) const {
		return &data_[limits_[i]];
	}

	const _t* end(size_t i) const {
		return &data_[limits_[i + 1]];
	}

	_t* begin(size_t i) {
		return &data_[limits_[i]];
	}

	_t* end(size_t i) {
		return &data_[limits_[i + 1]];
	}

private:

	std::vector<_t> data_;
	std::vector<size_t> limits_;

};

#endif