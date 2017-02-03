/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
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

	String_set(Input_stream &file)
	{
		file.read(limits_);
		file.read(data_);
	}

	static void skip(Input_stream &file)
	{
		file.skip_vector<size_t>();
		file.skip_vector<_t>();
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
