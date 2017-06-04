/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef STRING_SET_H_
#define STRING_SET_H_

#include <vector>
#include <stddef.h>
#include "../util/binary_file.h"

using std::vector;

template<char _pchar = '\xff', size_t _padding = 1lu>
struct String_set
{

	typedef char _t;
	enum { PERIMETER_PADDING = 256 };

	String_set():
		data_ (PERIMETER_PADDING)
	{ limits_.push_back(PERIMETER_PADDING); }

	void finish_reserve()
	{
		data_.resize(raw_len() + PERIMETER_PADDING);
		for(unsigned i=0;i<PERIMETER_PADDING;++i) {
			data_[i] = _pchar;
			data_[raw_len()+i] = _pchar;
		}
	}

	void reserve(size_t n)
	{
		limits_.push_back(raw_len() + n + _padding);
	}

	void push_back(const vector<_t> &v)
	{
		limits_.push_back(raw_len() + v.size() + _padding);
		data_.insert(data_.end(), v.begin(), v.end());
		data_.insert(data_.end(), _padding, _pchar);
	}

	void fill(size_t n, _t v)
	{
		limits_.push_back(raw_len() + n + _padding);
		data_.insert(data_.end(), n, v);
		data_.insert(data_.end(), _padding, _pchar);
	}

	_t* ptr(size_t i)
	{ return &data_[limits_[i]]; }

	const _t* ptr(size_t i) const
	{ return &data_[limits_[i]]; }

	size_t check_idx(size_t i) const
	{
		if (limits_.size() < i + 2)
			throw std::runtime_error("Sequence set index out of bounds.");
		return i;
	}

	size_t length(size_t i) const
	{ return limits_[i+1] - limits_[i] - _padding; }

	size_t get_length() const
	{ return limits_.size() - 1; }

	void save(Output_stream &file) const
	{
		file.write(limits_);
		file.write(data_);
	}

	size_t raw_len() const
	{ return limits_.back(); }

	size_t letters() const
	{ return raw_len() - get_length() - PERIMETER_PADDING; }

	_t* data(ptrdiff_t p = 0)
	{ return &data_[p]; }

	const _t* data(ptrdiff_t p = 0) const
	{ return &data_[p]; }

	size_t position(const _t* p) const
	{ return p - data(); }

	size_t position(size_t i, size_t j) const
	{ return limits_[i] + j; }

	std::pair<size_t,size_t> local_position(size_t p) const
	{
		size_t i = std::upper_bound(limits_.begin(), limits_.end(), p) - limits_.begin() - 1;
		return std::pair<size_t,size_t> (i, p - limits_[i]);
	}

	sequence operator[](size_t i) const
	{ return sequence (ptr(i), length(i)); }

private:

	vector<_t> data_;
	vector<size_t> limits_;

};

#endif /* STRING_SET_H_ */
