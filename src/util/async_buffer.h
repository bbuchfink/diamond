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

#ifndef ASYNC_BUFFER_H_
#define ASYNC_BUFFER_H_

#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>

namespace io = boost::iostreams;

using std::vector;
using std::string;
using std::endl;
using boost::ptr_vector;

struct Buffer_file_read_exception : public diamond_exception
{
	Buffer_file_read_exception(const char* file_name, size_t count, size_t n):
		diamond_exception (string("Error reading buffer file ") + file_name + " (" + boost::lexical_cast<string>(count) + '/' + boost::lexical_cast<string>(n) + ')')
	{ }
};

const unsigned async_buffer_max_bins = 4;

template<typename _t>
struct Async_buffer
{

	typedef vector<_t> Vector;

	Async_buffer(size_t input_count, const string &tmpdir, unsigned bins):
		bins_ (bins),
		bin_size_ ((input_count + bins_ - 1) / bins_)
	{
		log_stream << "Async_buffer() " << input_count << ',' << bin_size_ << endl;
		for(unsigned j=0;j<program_options::threads();++j)
			for(unsigned i=0;i<bins;++i) {
				tmp_file_.push_back(Temp_file ());
				out_.push_back(new Output_stream (tmp_file_.back()));
				size_.push_back(0);
			}
	}

	struct Iterator
	{
		Iterator(Async_buffer &parent, unsigned thread_num):
			parent_ (parent),
			thread_num_ (thread_num)
		{
			for(unsigned i=0;i<parent.bins_;++i) {
				size_[i] = 0;
				out_[i] = parent.get_out(thread_num_, i);
			}
		}
		void push(const _t &x)
		{
			const unsigned bin = x / parent_.bin_size_;
			assert(bin < parent_.bins());
			buffer_[bin*buffer_size+(size_[bin]++)] = x;
			if(size_[bin] == buffer_size)
				flush(bin);
		}
		void flush(unsigned bin)
		{
			out_[bin]->write(&buffer_[bin*buffer_size], size_[bin]);
			parent_.add_size(thread_num_, bin, size_[bin]);
			size_[bin] = 0;
		}
		~Iterator()
		{
			for(unsigned bin=0;bin<parent_.bins_;++bin)
				flush(bin);
		}
	private:
		enum { buffer_size = 65536 };
		_t buffer_[async_buffer_max_bins*buffer_size];
		size_t size_[async_buffer_max_bins];
		Output_stream* out_[async_buffer_max_bins];
		Async_buffer &parent_;
		const unsigned thread_num_;
	};

	void close()
	{
		for(ptr_vector<Output_stream>::iterator i=out_.begin();i!=out_.end();++i)
			i->close();
		out_.clear();
		log_stream << "Async_buffer.close() " << endl;
	}

	void load(vector<_t> &data, unsigned bin) const
	{
		size_t size = 0;
		for(unsigned i=0;i<program_options::threads();++i)
			size += size_[i*bins_+bin];
		log_stream << "Async_buffer.load() " << size << endl;
		data.resize(size);
		_t* ptr = data.data();
		for(unsigned i=0;i<program_options::threads();++i) {
			Input_stream f (tmp_file_[i*bins_+bin]);
			const size_t s = size_[i*bins_+bin];
			const size_t n = f.read(ptr, s);
			ptr += s;
			f.close();
			if(n != s)
				throw Buffer_file_read_exception(f.file_name.c_str(), s, n);
		}
	}

	unsigned bins() const
	{ return bins_; }

private:

	Output_stream* get_out(unsigned threadid, unsigned bin)
	{ return &out_[threadid*bins_+bin]; }

	void add_size(unsigned thread_id, unsigned bin, size_t n)
	{ size_[thread_id*bins_+bin] += n; }

	const unsigned bins_, bin_size_;
	ptr_vector<Output_stream> out_;
	vector<size_t> size_;
	vector<Temp_file> tmp_file_;

};

#endif /* ASYNC_BUFFER_H_ */
