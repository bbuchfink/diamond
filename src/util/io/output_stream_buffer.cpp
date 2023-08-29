/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include "../../basic/config.h"
#include "output_stream_buffer.h"

using std::pair;

OutputStreamBuffer::OutputStreamBuffer(StreamEntity* prev):
	StreamEntity(prev),
	buf_size_(prev->file_name().empty() ? STDOUT_BUF_SIZE : config.file_buffer_size),
	buf_(new char[buf_size_])
{}

pair<char*, char*> OutputStreamBuffer::write_buffer()
{
	return std::make_pair(buf_.get(), buf_.get() + buf_size_);
}

void OutputStreamBuffer::flush(size_t count)
{
	prev_->write(buf_.get(), count);
}

void OutputStreamBuffer::seek(int64_t p, int origin)
{
	prev_->seek(p, origin);
}

void OutputStreamBuffer::rewind()
{
	prev_->rewind();
}

int64_t OutputStreamBuffer::tell()
{
	return prev_->tell();
}
