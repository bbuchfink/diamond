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

#ifndef TRACE_PT_BUFFER_H_
#define TRACE_PT_BUFFER_H_

#include "../util/async_buffer.h"
#include "../basic/match.h"

using std::auto_ptr;

template<typename _locr, typename _locl>
struct Trace_pt_buffer : public Async_buffer<hit<_locr,_locl> >
{
	Trace_pt_buffer(size_t input_size, const string &tmpdir, bool mem_buffered):
		Async_buffer<hit<_locr,_locl> > (input_size, tmpdir, mem_buffered ? mem_bins : file_bins)
	{ }
	enum { mem_bins = 1, file_bins = 4 };
	static Trace_pt_buffer *instance;
};

template<typename _locr, typename _locl> Trace_pt_buffer<_locr,_locl>* Trace_pt_buffer<_locr,_locl>::instance;

template<typename _locr, typename _locl>
struct Trace_pt_list : public vector<hit<_locr,_locl> >
{
	void init()
	{
		pos_ = this->begin();
	}
	struct Query_range
	{
		Query_range(Trace_pt_list &parent):
			parent_ (parent)
		{ }
		bool operator()()
		{
			begin = parent_.pos_;
			end = std::min(begin + 4096, parent_.end());
			if(end >= parent_.end())
				return false;
			const unsigned c = query_contexts(), q = end->query_/c;
			for(; end<parent_.end() && end->query_/c == q; ++end);
			parent_.pos_ = end;
			return end < parent_.end();
		}
		typename Trace_pt_list::iterator begin, end;
	private:
		Trace_pt_list &parent_;
	};
	Query_range get_range()
	{ return Query_range (*this); }
private:
	typename vector<hit<_locr,_locl> >::iterator pos_;
};

#endif /* TRACE_PT_BUFFER_H_ */
