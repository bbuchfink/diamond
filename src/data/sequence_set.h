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

#ifndef SEQUENCE_SET_H_
#define SEQUENCE_SET_H_

#include <iostream>
#include <string>
#include "../basic/sequence.h"
#include "string_set.h"

using std::cout;
using std::endl;
using std::pair;

template<typename _val>
struct Sequence_set : public String_set<_val>
{

	Sequence_set()
	{ }

	Sequence_set(Input_stream &file):
		String_set<_val> (file)
	{ }

	void print_stats() const
	{ verbose_stream << "Sequences = " << this->get_length() << ", letters = " << this->letters() << endl; }

	pair<size_t,size_t> len_bounds(size_t min_len) const
	{
		const size_t l (this->get_length());
		size_t max = 0, min = std::numeric_limits<size_t>::max();
		for(size_t i=0;i<l;++i) {
			max = std::max(this->length(i), max);
			min = this->length(i) >= min_len ? std::min(this->length(i), min) : min;
		}
		return pair<size_t,size_t> (min, max);
	}

	sequence<const _val> window_infix(size_t offset, unsigned &left) const
	{
		const _val* begin (this->data(offset));
		unsigned n (0);
		while(*begin != String_set<_val>::PADDING_CHAR && n <= program_options::window) {
			--begin;
			++n;
		}
		++begin;
		left = program_options::window + 1 - n;
		const _val* end (this->data(offset));
		n = 0;
		while(*end != String_set<_val>::PADDING_CHAR && n < program_options::window) {
			++end;
			++n;
		}
		return sequence<const _val> (begin, end - begin);
	}

	sequence<const _val> fixed_window_infix(size_t offset) const
	{
		const _val* begin (this->data(offset));
		unsigned n (0);
		while(*begin != String_set<_val>::PADDING_CHAR && n <= program_options::window) {
			--begin;
			++n;
		}
		++begin;
		const _val* s (this->data(offset - program_options::window));
		return sequence<const _val> (s, 2*program_options::window, begin - s);
	}

	vector<size_t> partition() const
	{
		vector<size_t> v;
		const size_t l = (this->letters()+Const::seqp-1) / Const::seqp;
		v.push_back(0);
		for(unsigned i=0;i<this->get_length();) {
			size_t n = 0;
			while(i<this->get_length() && n < l)
				n += this->length(i++);
			v.push_back(i);
		}
		for(unsigned i=v.size();i<Const::seqp+1;++i)
			v.push_back(this->get_length());
		return v;
	}

	size_t reverse_translated_len(size_t i) const
	{
		const size_t j (i - i%6);
		const size_t l (this->length(j));
		if(this->length(j+2) == l)
			return l*3 + 2;
		else if(this->length(j+1) == l)
			return l*3 + 1;
		else
			return l*3;
	}

	virtual ~Sequence_set()
	{ }

};

#endif /* SEQUENCE_SET_H_ */
