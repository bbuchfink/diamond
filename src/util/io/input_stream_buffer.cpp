/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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

#include <utility>
#include "input_stream_buffer.h"
#include "basic/config.h"

using std::make_pair;
using std::pair;

InputStreamBuffer::InputStreamBuffer(StreamEntity* prev, int flags) :
	StreamEntity(prev, prev->seekable()),
	begin(nullptr),
	end(nullptr),
	buf_size_(config.file_buffer_size),
	buf_(new char[buf_size_]),
	load_buf_((flags& ASYNC) != 0 ? new char[buf_size_] : nullptr),
	load_count_(0),
	async_((flags & ASYNC) != 0),
	load_worker_(nullptr)
{
}

void InputStreamBuffer::rewind()
{
	prev_->rewind();
	file_offset_ = 0;
	begin = end = nullptr;
}

void InputStreamBuffer::seek(int64_t pos, int origin)
{
	prev_->seek(pos, origin);
	file_offset_ = 0;
	begin = end = nullptr;
}

bool InputStreamBuffer::fetch()
{
	if (load_worker_) {
		load_worker_->join();
		delete load_worker_;
		load_worker_ = nullptr;
		std::swap(buf_, load_buf_);
		begin = buf_.get();
		end = begin + load_count_;
	}
	else {
		const size_t n = prev_->read(buf_.get(), buf_size_);
		if (prev_->seekable())
			file_offset_ = prev_->tell();
		begin = buf_.get();
		end = begin + n;
	}
	if (async_)
		load_worker_ = new std::thread(load_worker, this);	
	return begin != nullptr && end > begin;
}

void InputStreamBuffer::load_worker(InputStreamBuffer* buf) {
	buf->load_count_ = buf->prev_->read(buf->load_buf_.get(), buf->buf_size_);
}

void InputStreamBuffer::close() {
	if (load_worker_) {
		load_worker_->join();
		delete load_worker_;
		load_worker_ = nullptr;
	}
	prev_->close();
}

int64_t InputStreamBuffer::tell() {
	if (!seekable())
		throw std::runtime_error("Calling tell on non seekable stream.");
	return file_offset_;
}

bool InputStreamBuffer::eof() {
	return prev_->eof();
}