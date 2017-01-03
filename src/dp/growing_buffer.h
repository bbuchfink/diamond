/****
Copyright (c) 2014, University of Tuebingen
Author: Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#ifndef GROWING_BUFFER_H_
#define GROWING_BUFFER_H_

#include <vector>

using std::vector;
using std::pair;

template<typename _t>
struct Growing_buffer
{

	inline void init(size_t size, size_t padding, size_t padding_front, _t init)
	{
		const size_t total = size + padding;
		data_.clear();
		data_.resize(total);
		col_size_ = total;
		for(size_t i=0;i<total;++i)
			data_[i] = init;
		init_ = init;
		center_.clear();
		center_.push_back(-1);
	}

	inline pair<_t*,_t*> get(int center)
	{
		data_.resize(data_.size() + col_size_);
		_t* ptr = last();
		for(size_t i=0;i<col_size_;++i)
			ptr[i] = init_;
		center_.push_back(center);
		return pair<_t*,_t*> (ptr-col_size_, ptr);
	}

	inline _t* last()
	{ return &*(data_.end() - col_size_); }

	const _t* column(int col) const
	{ return &data_[col_size_*col]; }

	int center(int col) const
	{ return center_[col]; }

private:
	vector<_t> data_;
	vector<int> center_;
	size_t col_size_;
	_t init_;

};

#endif
