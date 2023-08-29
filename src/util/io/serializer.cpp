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

#include <algorithm>
#include <string.h>
#include "serializer.h"

using std::string;
using std::pair;

Serializer::Serializer(StreamEntity *buffer, int flags) :
	buffer_(buffer),
	varint_(flags & VARINT)
{
	reset_buffer();
}

void Serializer::flush()
{
	buffer_->flush(next_ - begin_);
}

void Serializer::close()
{
	flush();
	buffer_->close();
}

void Serializer::seek(int64_t p, int origin)
{
	flush();
	buffer_->seek(p, origin);
	reset_buffer();
}

void Serializer::rewind()
{
	flush();
	buffer_->rewind();
	reset_buffer();
}

size_t Serializer::tell()
{
	flush();
	reset_buffer();
	return buffer_->tell();
}

string Serializer::file_name() const
{
	return buffer_->file_name();
}

FILE* Serializer::file()
{
	return buffer_->file();
}

Serializer::~Serializer()
{
	delete buffer_;
}

void Serializer::reset_buffer()
{
	pair<char*, char*> buf = buffer_->write_buffer();
	begin_ = next_ = buf.first;
	end_ = buf.second;
}

void Serializer::write_raw(const char *ptr, size_t count)
{
	do {
		size_t n = std::min(avail(), count);
		memcpy(next_, ptr, n);
		next_ += n;
		ptr += n;
		count -= n;
		if (avail() == 0) {
			flush();
			reset_buffer();
		}
	} while (count > 0);
}

void Serializer::set(int flag)
{
	switch (flag) {
	case VARINT:
		varint_ = true;
	default:
		;
	}
}

void Serializer::unset(int flag)
{
	switch (flag) {
	case VARINT:
		varint_ = false;
	default:
		;
	}
}

void Serializer::consume(const char *ptr, size_t n) {
	write_raw(ptr, n);
}

void Serializer::finalize() {
	close();
}