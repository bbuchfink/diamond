/****
Copyright (c) 2015, University of Tuebingen
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

#ifndef PACKED_SEQUENCE_H_
#define PACKED_SEQUENCE_H_

#include <vector>
#include "../util/binary_buffer.h"

using std::vector;

bool has_n(const sequence<const Nucleotide> &seq)
{
	for(unsigned i=0;i<seq.length();++i)
		if(seq[i] == Value_traits<Nucleotide>::MASK_CHAR)
			return true;
	return false;
}

struct Packed_sequence
{

	Packed_sequence(const sequence<const Nucleotide> &seq):
		has_n_ (::has_n(seq))
	{
		if(has_n_)
			pack<Nucleotide,3>(seq);
		else
			pack<Nucleotide,2>(seq);
	}

	Packed_sequence(const sequence<const Amino_acid> &seq):
		has_n_ (false)
	{ pack<Amino_acid,5>(seq); }

	Packed_sequence(Binary_buffer::Iterator &it, unsigned len, bool has_n, unsigned b):
		has_n_ (has_n)
	{
		const size_t l = (len*b+7)/8;
		it.read(data_, l);
	}

	template<typename _val>
	void unpack(vector<_val> &dst, unsigned b, unsigned len)
	{
		dst.clear();
		unsigned x = 0, n = 0, l = 0;
		const unsigned mask = (1<<b)-1;
		for(unsigned i=0;i<data_.size();++i) {
			x |= (unsigned)data_[i] << n;
			n += 8;
			while(n >= b && l < len) {
				dst.push_back(x & mask);
				n -= b;
				x >>= b;
				++l;
			}
		}
	}

	const vector<uint8_t>& data() const
	{ return data_; }

	bool has_n() const
	{ return has_n_; }

private:

	template<typename _val, unsigned _b>
	void pack(const sequence<const _val> &seq)
	{
		unsigned x = 0, n = 0;
		for(unsigned i=0;i<seq.length();++i) {
			x |= (unsigned)seq[i] << n;
			n += _b;
			if(n >= 8) {
				data_.push_back(x & 0xff);
				n -= 8;
				x >>= 8;
			}
		}
		if(n > 0)
			data_.push_back(x & 0xff);
	}

	bool has_n_;
	vector<uint8_t> data_;

};

#endif /* PACKED_SEQUENCE_H_ */
