/****
Copyright (c) 2013-2017, Benjamin Buchfink, University of Tuebingen
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

#ifndef ASYNC_BUFFER_H_
#define ASYNC_BUFFER_H_

#include <vector>
#include <exception>
#include "../basic/config.h"
#include "temp_file.h"

using std::vector;
using std::string;
using std::endl;

template<typename _t>
struct Async_buffer
{

	typedef vector<_t> Vector;

	Async_buffer(size_t input_count, const string &tmpdir, unsigned bins) :
		bins_(bins),
		bin_size_((input_count + bins_ - 1) / bins_)
	{
		log_stream << "Async_buffer() " << input_count << ',' << bin_size_ << endl;
		for (unsigned j = 0; j < config.threads_; ++j)
			for (unsigned i = 0; i < bins; ++i) {
				tmp_file_.push_back(Temp_file());
				size_.push_back(0);
			}
	}

	struct Iterator
	{
		Iterator(Async_buffer &parent, unsigned thread_num) :
			buffer_(parent.bins()),
			parent_(parent),
			thread_num_(thread_num)
		{
			for (unsigned i = 0; i < parent.bins_; ++i)
				out_.push_back(parent.get_out(thread_num_, i));
		}
		void push(const _t &x)
		{
			const unsigned bin = (unsigned)(x / parent_.bin_size_);
			assert(bin < parent_.bins());
			buffer_[bin].push_back(x);
			if (buffer_[bin].size() == buffer_size)
				flush(bin);
		}
		void flush(unsigned bin)
		{
			out_[bin]->typed_write(buffer_[bin].data(), buffer_[bin].size());
			parent_.add_size(thread_num_, bin, buffer_[bin].size());
			buffer_[bin].clear();
		}
		~Iterator()
		{
			for (unsigned bin = 0; bin < parent_.bins_; ++bin)
				flush(bin);
		}
	private:
		enum { buffer_size = 65536 };
		vector<vector<_t> > buffer_;
		vector<Temp_file*> out_;
		Async_buffer &parent_;
		const unsigned thread_num_;
	};

	size_t load(vector<_t> &data, unsigned bin) const
	{
		static size_t total_size;
		if (bin == 0)
			total_size = 0;
		size_t size = 0;
		for (unsigned i = 0; i < config.threads_; ++i)
			size += size_[i*bins_ + bin];
		log_stream << "Async_buffer.load() " << size << "(" << (double)size*sizeof(_t) / (1 << 30) << " GB)" << endl;
		total_size += size;
		data.resize(size);
		_t* ptr = data.data();
		for (unsigned i = 0; i < config.threads_; ++i) {
			Input_stream f(tmp_file_[i*bins_ + bin]);
			const size_t s = size_[i*bins_ + bin];
			const size_t n = f.read(ptr, s);
			ptr += s;
			f.close_and_delete();
			if (n != s)
				throw std::runtime_error("Error reading temporary file: " + f.file_name);
		}
		return total_size*sizeof(_t);
	}

	unsigned bins() const
	{
		return bins_;
	}

private:

	Temp_file* get_out(unsigned threadid, unsigned bin)
	{
		return &tmp_file_[threadid*bins_ + bin];
	}

	void add_size(unsigned thread_id, unsigned bin, size_t n)
	{
		size_[thread_id*bins_ + bin] += n;
	}

	const unsigned bins_;
	const size_t bin_size_;
	vector<size_t> size_;
	vector<Temp_file> tmp_file_;

};

#endif /* ASYNC_BUFFER_H_ */
