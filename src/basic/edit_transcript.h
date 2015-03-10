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

#ifndef EDIT_TRANSCRIPT_H_
#define EDIT_TRANSCRIPT_H_

#include <vector>
#include "packed_transcript.h"

using std::vector;
using std::endl;

struct Edit_transcript : public vector<char>
{

	void push(Edit_operation op, unsigned n)
	{
		for(unsigned i=0;i<n;++i)
			this->push_back(op);
	}

	template<typename _val>
	void print_packed(Text_buffer& buf, const sequence<const _val> &query, const sequence<const _val> &subject, unsigned qpos, unsigned spos)
	{
		for(const_iterator i = this->end()-1; i >= this->begin();)
			switch(*i) {
			case op_match:
				if(query[qpos] == mask_critical(subject[spos]))
					print_match(buf, i, query, subject, qpos, spos);
				else {
					buf.write(Packed_operation(op_substitution, mask_critical(subject[spos])));
					++qpos;
					++spos;
					--i;
				}
				break;
			case op_insertion:
				print_insertion(buf, i, qpos);
				break;
			case op_deletion:
				print_deletion(buf, i, subject, spos);
			}
		buf.write(Packed_operation::terminator());
	}

	template<typename _val>
	void print_match(Text_buffer& buf, const_iterator& i, const sequence<const _val> &query, const sequence<const _val>& subject, unsigned& qpos, unsigned &spos)
	{
		unsigned n=0;
		for(;i >= this->begin() && *i == op_match && query[qpos] == mask_critical(subject[spos]); --i) {
			++qpos;
			++spos;
			++n;
		}
		print_number(buf, n, op_match);
	}

	template<typename _val>
	void print_deletion(Text_buffer& buf, const_iterator& i, const sequence<const _val>& subject, unsigned &spos)
	{
		for(;i >= this->begin() && *i == op_deletion; --i)
			buf.write(Packed_operation(op_deletion, mask_critical(subject[spos++])));
	}

	void print_insertion(Text_buffer &buf, const_iterator& i, unsigned &qpos)
	{
		unsigned n = 0;
		for(;i >= this->begin() && *i == op_insertion; --i) {
			++n;
			++qpos;
		}
		print_number(buf, n, op_insertion);
	}

	void print_number(Text_buffer &buf, unsigned n, Edit_operation op)
	{
		while(n>0) {
			unsigned m = std::min(n, 63u);
			buf.write(Packed_operation(op, m));
			n -= m;
		}
	}

	template<typename _val>
	void print(std::ostream &os, const _val *query, const _val *subject)
	{
		print(os, query, op_deletion);
		os << endl;
		print(os, subject, op_insertion);
		os << endl;
	}

	template<typename _val>
	void print(std::ostream &os, const _val *s, Edit_operation gap_op)
	{
		for(const_iterator i = this->end()-1; i >= this->begin(); --i)
			if(*i == gap_op)
				os << '-';
			else
				os << Value_traits<_val>::ALPHABET[*(s++)];
	}

private:

	void print_matches(char *&ptr, unsigned &n)
	{
		if(n > 0) {
			ptr += sprintf(ptr, "%u", n);
			n = 0;
		}
	}

};

#endif /* EDIT_TRANSCRIPT_H_ */
