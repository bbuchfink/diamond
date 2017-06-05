/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef ASYNC_BUFFER_H_
#define ASYNC_BUFFER_H_

#include <vector>
#include <exception>
#include "../basic/config.h"
#include "temp_file.h"
#include "log_stream.h"

using std::vector;
using std::string;
using std::endl;

template<typename _t>
struct Async_buffer
{

	typedef vector<_t> Vector;

	Async_buffer(size_t input_count, const string &tmpdir, unsigned bins) :
		bins_(bins),
		bin_size_((input_count + bins_ - 1) / bins_),
		input_count_(input_count),
		bins_processed_(0)
	{
		log_stream << "Async_buffer() " << input_count << ',' << bin_size_ << endl;
		for (unsigned j = 0; j < config.threads_; ++j)
			for (unsigned i = 0; i < bins; ++i) {
				tmp_file_.push_back(Temp_file());
				size_.push_back(0);
			}
	}

	size_t begin(size_t bin) const
	{
		return bin*bin_size_;
	}

	size_t end(size_t bin) const
	{
		return std::min((bin + 1)*bin_size_, input_count_);
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
			out_[bin]->write(buffer_[bin].data(), buffer_[bin].size());
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

	size_t load(vector<_t> &data, size_t max_size, std::pair<size_t,size_t> &input_range)
	{
		static size_t total_size;
		if (bins_processed_ == 0)
			total_size = 0;
		if (bins_processed_ == bins_) {
			input_range = std::make_pair(0, 0);
			return total_size*sizeof(_t);
		}
		size_t size = bin_size(bins_processed_), end = bins_processed_ + 1, current_size;
		while (end < bins_ && (size + (current_size = bin_size(end)))*sizeof(_t) < max_size) {
			size += current_size;
			++end;
		}
		log_stream << "Async_buffer.load() " << size << "(" << (double)size*sizeof(_t) / (1 << 30) << " GB)" << endl;
		total_size += size;
		data.resize(size);
		_t* ptr = data.data();
		input_range.first = begin(bins_processed_);
		for (; bins_processed_ < end; ++bins_processed_)
			load_bin(ptr, bins_processed_);
		input_range.second = this->end(bins_processed_ - 1);
		return total_size*sizeof(_t);
	}

	unsigned bins() const
	{
		return bins_;
	}

private:

	void load_bin(_t*& ptr, size_t bin)
	{
		for (unsigned i = 0; i < config.threads_; ++i) {
			Input_stream f(tmp_file_[i*bins_ + bin]);
			const size_t s = size_[i*bins_ + bin];
			const size_t n = f.read(ptr, s);
			ptr += s;
			f.close_and_delete();
			if (n != s)
				throw std::runtime_error("Error reading temporary file: " + f.file_name);
		}
	}

	size_t bin_size(size_t bin) const
	{
		size_t size = 0;
		for (unsigned i = 0; i < config.threads_; ++i)
			size += size_[i*bins_ + bin];
		return size;
	}

	Temp_file* get_out(unsigned threadid, unsigned bin)
	{
		return &tmp_file_[threadid*bins_ + bin];
	}

	void add_size(unsigned thread_id, unsigned bin, size_t n)
	{
		size_[thread_id*bins_ + bin] += n;
	}

	const unsigned bins_;
	const size_t bin_size_, input_count_;
	size_t bins_processed_;
	vector<size_t> size_;
	vector<Temp_file> tmp_file_;

};

#endif /* ASYNC_BUFFER_H_ */
