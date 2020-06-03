/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include <exception>
#include <assert.h>
#include "../basic/config.h"
#include "io/temp_file.h"
#include "io/input_file.h"
#include "log_stream.h"
#include "../util/ptr_vector.h"
#include "io/async_file.h"

template<typename _t>
struct Async_buffer
{

	typedef std::vector<_t> Vector;

	Async_buffer(size_t input_count, const std::string &tmpdir, unsigned bins) :
		bins_(bins),
		bin_size_((input_count + bins_ - 1) / bins_),
		input_count_(input_count),
		bins_processed_(0)
	{
		log_stream << "Async_buffer() " << input_count << ',' << bin_size_ << endl;
		for (unsigned i = 0; i < bins; ++i)
			tmp_file_.push_back(new AsyncFile());
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
		Iterator(Async_buffer &parent, size_t thread_num) :
			buffer_(parent.bins()),
			parent_(parent)
		{
			for (unsigned i = 0; i < parent.bins_; ++i)
				out_.push_back(&parent.tmp_file_[i]);
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
			buffer_[bin].clear();
		}
		~Iterator()
		{
			for (unsigned bin = 0; bin < parent_.bins_; ++bin)
				flush(bin);
		}
	private:
		enum { buffer_size = 65536 };
		std::vector<std::vector<_t>> buffer_;
		std::vector<AsyncFile*> out_;
		Async_buffer &parent_;
	};

	size_t load(std::vector<_t> &data, size_t max_size, std::pair<size_t,size_t> &input_range)
	{
		static size_t total_size;
		if (bins_processed_ == 0)
			total_size = 0;
		if (bins_processed_ == bins_) {
			input_range = std::make_pair(0, 0);
			return total_size*sizeof(_t);
		}
		size_t size = tmp_file_[bins_processed_].tell()/sizeof(_t), end = bins_processed_ + 1, current_size;
		while (end < bins_ && (size + (current_size = tmp_file_[end].tell() / sizeof(_t))) * sizeof(_t) < max_size) {
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
		const size_t s = tmp_file_[bin].tell() / sizeof(_t);
		InputFile f(tmp_file_[bin]);
		const size_t n = f.read(ptr, s);
		ptr += s;
		f.close_and_delete();
		if (n != s)
			throw std::runtime_error("Error reading temporary file: " + f.file_name);
	}

	const unsigned bins_;
	const size_t bin_size_, input_count_;
	size_t bins_processed_;
	PtrVector<AsyncFile> tmp_file_;

};