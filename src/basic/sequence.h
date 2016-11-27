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

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <iostream>
#include <vector>
#include "../basic/value.h"
#include "../util/binary_buffer.h"
#include "../util/text_buffer.h"

using std::vector;

struct sequence
{
	struct Reversed {};
	sequence():
		len_ (0),
		clipping_offset_ (0),
		data_ (0)
	{ }
	sequence(const Letter *data, size_t len, int clipping_offset = 0):
		len_ (len),
		clipping_offset_ (clipping_offset),
		data_ (data)
	{ }
	sequence(const vector<Letter> &data):
		len_(data.size()),
		clipping_offset_(0),
		data_(data.data())
	{}
	sequence(const sequence &seq, int from, int to):
		len_(to-from+1),
		clipping_offset_(0),
		data_(&seq[from])
	{}
	size_t length() const
	{
		return len_;
	}
	size_t clipped_length() const
	{
		return len_ - clipping_offset_;
	}
	size_t aligned_clip(unsigned padding) const
	{
		return clipping_offset_ > (int)padding ? clipping_offset_ - padding : 0;
	}
	const Letter* data() const
	{
		return data_;
	}
	const Letter* clipped_data() const
	{
		return data_ + clipping_offset_;
	}
	const Letter* aligned_data(unsigned padding) const
	{
		return data_ + padding;
	}
	const Letter& operator [](size_t i) const
	{
		return data_[i];
	}
	bool empty() const
	{ return len_ == 0; }
	const char* c_str() const
	{ return reinterpret_cast<const char*>(data_); }
	size_t print(char *ptr, unsigned begin, unsigned len) const
	{
		for(unsigned i=begin;i<begin+len;++i)
			*(ptr++) = to_char(data_[i]);
		return len;
	}
	Text_buffer& print(Text_buffer &buf, size_t begin, size_t end, const Value_traits& value_traits) const
	{
		for (size_t i = begin; i < end; ++i)
			buf << value_traits.alphabet[(long)data_[i]];
		return buf;
	}
	std::ostream& print(std::ostream &os, const Value_traits &v) const
	{
		for (unsigned i = 0; i<len_; ++i)
			os << v.alphabet[(long)data_[i]];
		return os;
	}
	std::ostream& print(std::ostream &os, const Value_traits &v, Reversed) const
	{
		for (int i = (int)len_ - 1; i >= 0; --i)
			os << v.alphabet[(long)data_[i]];
		return os;
	}
	sequence subseq(int begin, int end) const
	{
		return sequence(*this, begin, end - 1);
	}
	friend std::ostream& operator<<(std::ostream &os, const sequence &s)
	{
		return s.print(os, value_traits);
	}
	friend Text_buffer& operator<<(Text_buffer &buf, const sequence &s)
	{
		for(unsigned i=0;i<s.len_;++i)
			buf << value_traits.alphabet[(long)s.data_[i]];
		return buf;
	}
	/*friend std::ostream& operator<<(std::ostream &os, const sequence &s)
	{
		std::cout << "co = " << s.clipping_offset_ << std::endl;
		for(unsigned i=s.clipping_offset_;i<s.len_;++i) {
			if(s.data_[i] == 24)
				break;
			os << mask_critical(s.data_[i]);
		}
		return os;
	}*/
	static vector<Letter> from_string(const char* str);
	size_t			len_;
	int				clipping_offset_;
	const Letter	*data_;
};


#endif /* SEQUENCE_H_ */
