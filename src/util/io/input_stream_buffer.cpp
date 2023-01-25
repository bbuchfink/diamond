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

#include <algorithm>
#include <string.h>
#include <utility>
#include "input_stream_buffer.h"
#include "../../basic/config.h"

using std::make_pair;
using std::pair;

InputStreamBuffer::InputStreamBuffer(StreamEntity* prev, int flags) :
	StreamEntity(prev, prev->seekable()),
	buf_size_(config.file_buffer_size),
	buf_(new char[buf_size_]),
	load_buf_((flags& ASYNC) != 0 ? new char[buf_size_] : nullptr),
	putback_count_(0),
	load_count_(0),
	async_((flags & ASYNC) != 0),
	load_worker_(nullptr)
{
}

void InputStreamBuffer::rewind()
{
	prev_->rewind();
	file_offset_ = 0;
	putback_count_ = 0;
}

void InputStreamBuffer::seek(int64_t pos, int origin)
{
	prev_->seek(pos, origin);
	file_offset_ = 0;
}

void InputStreamBuffer::seek_forward(size_t n)
{
	prev_->seek_forward(n);
	file_offset_ = 0;
}

pair<const char*, const char*> InputStreamBuffer::read()
{
	size_t n;
	if (putback_count_ > 0) {
		n = putback_count_;
		putback_count_ = 0;
	}
	else {
		if (load_worker_) {
			load_worker_->join();
			delete load_worker_;
			load_worker_ = nullptr;
			std::swap(buf_, load_buf_);
			n = load_count_;
		}
		else {
			n = prev_->read(buf_.get(), buf_size_);
			if (prev_->seekable())
				file_offset_ = prev_->tell();
		}
		if (async_)
			load_worker_ = new std::thread(load_worker, this);
	}
	return make_pair(buf_.get(), buf_.get() + n);
}

void InputStreamBuffer::load_worker(InputStreamBuffer* buf) {
	buf->load_count_ = buf->prev_->read(buf->load_buf_.get(), buf->buf_size_);
}

void InputStreamBuffer::putback(const char* p, size_t n) {
	std::copy(p, p + n, buf_.get());
	putback_count_ = n;
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