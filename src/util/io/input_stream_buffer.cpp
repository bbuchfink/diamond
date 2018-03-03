/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include <algorithm>
#include <string.h>
#include "input_stream_buffer.h"

using std::make_pair;

InputStreamBuffer::InputStreamBuffer(StreamEntity *prev):
	StreamEntity(prev)
{
}

void InputStreamBuffer::rewind()
{
	prev_->rewind();
}

void InputStreamBuffer::seek(size_t pos)
{
	prev_->seek(pos);
}

void InputStreamBuffer::seek_forward(size_t n)
{
	prev_->seek_forward(n);
}

pair<const char*, const char*> InputStreamBuffer::read()
{
	size_t n = prev_->read(buf_, BUF_SIZE);
	return make_pair(buf_, buf_ + n);
}