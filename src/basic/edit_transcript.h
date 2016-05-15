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
#include "../util/util.h"

using std::vector;
using std::endl;

struct Edit_transcript
{

	Edit_transcript():
		begin_ (0),
		end_ (0)
	{ }

	Edit_transcript(const vector<char> &buf):
		begin_ (buf.size()),
		end_ (begin_)
	{ }

	Edit_transcript& set_end(const vector<char> &buf)
	{
		end_ = buf.size();
		return *this;
	}

	vector<char>::const_iterator begin(const vector<char> &buf) const
	{ return buf.begin() + begin_; }

	vector<char>::const_iterator end(const vector<char> &buf) const
	{ return buf.begin() + end_; }

private:

	void print_matches(char *&ptr, unsigned &n)
	{
		if(n > 0) {
			ptr += sprintf(ptr, "%u", n);
			n = 0;
		}
	}

	size_t begin_, end_;

};

struct Link_iterator
{
	Link_iterator(const Edit_transcript &right,
		const Edit_transcript &left,
		const vector<char> &transcript_buf):
			good_(true),
			i_ (left.begin(transcript_buf)),
			left_end_ (left.end(transcript_buf)),
			right_begin_ (right.begin(transcript_buf)),
			right_end_ (right.end(transcript_buf))
	{
		if(i_ >= left_end_-1)
			i_ = right_end_-1;
	}
	char operator*() const
	{ return *i_; }
	bool good() const
	//{ return i_ != right_begin_-1; }
	{
		return good_;
	}
	Link_iterator& operator++()
	{
		if (i_ == right_begin_)
			good_ = false;
		else if(i_ > right_begin_ && i_ < right_end_)
			--i_;
		else {
			++i_;
			if(i_ >= left_end_-1)
				i_ = right_end_-1;
		}
		return *this;
	}
private:
	bool good_;
	vector<char>::const_iterator i_;
	const vector<char>::const_iterator left_end_, right_begin_, right_end_;
};

inline void print_number(Text_buffer &buf, unsigned n, Edit_operation op)
{
	while(n>0) {
		unsigned m = std::min(n, 63u);
		buf.write(Packed_operation(op, m));
		n -= m;
	}
}

inline void print_match(Text_buffer& buf, Link_iterator& i, const sequence &query, const sequence& subject, unsigned& qpos, unsigned &spos)
{
	unsigned n=0;
	for(;i.good() && *i == op_match && query[qpos] == mask_critical(subject[spos]); ++i) {
		++qpos;
		++spos;
		++n;
	}
	print_number(buf, n, op_match);
}

inline void print_deletion(Text_buffer& buf, Link_iterator& i, const sequence& subject, unsigned &spos)
{
	for(;i.good() && *i == op_deletion; ++i)
		buf.write(Packed_operation(op_deletion, mask_critical(subject[spos++])));
}

inline void print_insertion(Text_buffer &buf, Link_iterator& i, unsigned &qpos)
{
	unsigned n = 0;
	for(;i.good() && *i == op_insertion; ++i) {
		++n;
		++qpos;
	}
	print_number(buf, n, op_insertion);
}

inline void print_packed(const Edit_transcript &right,
		const Edit_transcript &left,
		const vector<char> &transcript_buf,
		Text_buffer& buf,
		const sequence &query,
		const sequence &subject,
		unsigned qpos,
		unsigned spos)
{
	Link_iterator i (right, left, transcript_buf);
	while(i.good())
		switch(*i) {
		case op_match:
			if(query[qpos] == mask_critical(subject[spos]))
				print_match(buf, i, query, subject, qpos, spos);
			else {
				buf.write(Packed_operation(op_substitution, mask_critical(subject[spos])));
				++qpos;
				++spos;
				++i;
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

inline void print(Link_iterator i, std::ostream &os, const Letter *s, Edit_operation gap_op)
{
	unsigned n=0;
	for(; i.good(); ++i) {
		os << n << ' ';
		if(*i == gap_op)
			os << '-';
		else
			os << value_traits.alphabet[(long)*(s++)];
		++n;
	}
}

inline void print(std::ostream &os,
		const Letter *query,
		const Letter *subject,
		const Edit_transcript &right,
		const Edit_transcript &left,
		const vector<char> &transcript_buf)
{
	Link_iterator i (right, left, transcript_buf);
	print(i, os, query, op_deletion);
	os << endl;
	print(i, os, subject, op_insertion);
	os << endl;
}

#endif /* EDIT_TRANSCRIPT_H_ */
