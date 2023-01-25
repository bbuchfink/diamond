/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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
#include "value.h"
#include "../util/binary_buffer.h"
#include "sequence.h"

inline bool has_n(const Sequence &seq)
{
	for(Loc i=0;i<seq.length();++i)
		if(seq[i] == 4)
			return true;
	return false;
}

struct Packed_sequence
{

	Packed_sequence(const Sequence &seq, SequenceType type):
		has_n_ (type == SequenceType::nucleotide ? ::has_n(seq) : false)
	{
		switch (type) {
		case SequenceType::nucleotide:
			if (has_n_)
				pack<3>(seq);
			else
				pack<2>(seq);
			break;
		case SequenceType::amino_acid:
			pack<5>(seq);
		}
	}

	Packed_sequence(BinaryBuffer::Iterator &it, unsigned len, bool has_n, unsigned b):
		has_n_ (has_n)
	{
		const size_t l = (len*b+7)/8;
		it.read(data_, l);
	}

	void unpack(std::vector<Letter> &dst, unsigned b, unsigned len)
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

	const std::vector<uint8_t>& data() const
	{ return data_; }

	bool has_n() const
	{ return has_n_; }

private:

	template<unsigned _b>
	void pack(const Sequence &seq)
	{
		unsigned x = 0, n = 0;
		for(Loc i=0;i<seq.length();++i) {
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
	std::vector<uint8_t> data_;

};