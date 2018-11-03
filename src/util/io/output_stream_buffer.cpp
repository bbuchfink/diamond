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

#include "output_stream_buffer.h"

OutputStreamBuffer::OutputStreamBuffer(StreamEntity* prev):
	StreamEntity(prev)
{}

pair<char*, char*> OutputStreamBuffer::write_buffer()
{
	return std::make_pair(buf_, buf_ + BUF_SIZE);
}

void OutputStreamBuffer::flush(size_t count)
{
	prev_->write(buf_, count);
}

void OutputStreamBuffer::seek(size_t pos)
{
	prev_->seek(pos);
}

void OutputStreamBuffer::rewind()
{
	prev_->rewind();
}

size_t OutputStreamBuffer::tell()
{
	return prev_->tell();
}
