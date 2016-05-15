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

struct Trace_pt_buffer : public Async_buffer<hit>
{
	Trace_pt_buffer(size_t input_size, const string &tmpdir, bool mem_buffered):
		Async_buffer<hit> (input_size, tmpdir, mem_buffered ? mem_bins : file_bins)
	{ }
	enum { mem_bins = 1, file_bins = 4 };
	static Trace_pt_buffer *instance;
};

struct Trace_pt_list : public vector<hit>
{
	void init()
	{
		pos_ = this->begin();
		total_ = 0;
		count_ = 1;
#ifdef PRE_PARTITION
		p_.clear();
		p_.push_back(0);
		idx_ = 0;
		const unsigned c = query_contexts();
		typename vector<hit<_locr,_locl> >::iterator i = this->begin();
		unsigned total=0,count=1;
		for(; i < this->end();) {
			unsigned n=0;
			const unsigned min_size = std::max(4*total/count/5 + 1, program_options::fetch_size);
			for(;i<this->end() && n<min_size;) {
				const unsigned q = i->query_/c;
				for(; i<this->end() && i->query_/c == q; ++i)
					++n;
			}
			++count;
			total += n;
			p_.push_back(i - this->begin());
		}
		p_.push_back(i - this->begin());
#endif
	}
	struct Query_range
	{
		Query_range(Trace_pt_list &parent):
			parent_ (parent)
		{ }
#ifndef PRE_PARTITION
		bool operator()()
		{

			begin = parent_.pos_;
			//end = std::min(std::max(begin + 3*parent_.total_/parent_.count_/4 + 1, begin+program_options::fetch_size), parent_.end());
#ifdef NDEBUG
			end = std::min(begin + 3*parent_.total_/parent_.count_/4 + 1, parent_.end());
#else
			//end = parent_.end();
			ptrdiff_t x = std::min((ptrdiff_t)(3 * parent_.total_ / parent_.count_ / 4 + 1), parent_.end() - begin);
			end = std::min(begin + x, parent_.end());
#endif
			if(end >= parent_.end())
				return false;
			const unsigned c = query_contexts(), q = end->query_/c;
			for(; end<parent_.end() && end->query_/c == q; ++end);
			parent_.pos_ = end;
			parent_.total_ += end - begin;
			++parent_.count_;
			return end < parent_.end();
		}
#else
		bool operator()()
		{
			begin = parent_.begin()+parent_.p_[parent_.idx_];
			end = parent_.begin()+parent_.p_[parent_.idx_+1];
			printf("%lu %lu %lu\n", parent_.p_[parent_.idx_], parent_.p_[parent_.idx_+1], parent_.p_[parent_.idx_+1]-parent_.p_[parent_.idx_]);
			++parent_.idx_;
			return parent_.idx_ < parent_.p_.size()-1;
		}
#endif
		Trace_pt_list::iterator begin, end;
	private:
		Trace_pt_list &parent_;
	};
	Query_range get_range()
	{ return Query_range (*this); }
private:
	typename vector<hit>::iterator pos_;
#ifdef PRE_PARTITION
	vector<size_t> p_;
	unsigned idx_;
#else
	size_t total_, count_;
#endif
};

#endif /* TRACE_PT_BUFFER_H_ */

