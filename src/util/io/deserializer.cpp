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

#include <assert.h>
#include "deserializer.h"

using std::string;
using std::runtime_error;

Deserializer::Deserializer(InputStreamBuffer* buffer):
	buffer_(buffer)
{
}

Deserializer::~Deserializer()
{
	delete buffer_;
}

void Deserializer::close()
{
	buffer_->close();
}

void Deserializer::rewind()
{
	buffer_->rewind();
}

Deserializer& Deserializer::seek(int64_t pos)
{
	if (buffer_->seekable() && buffer_->tell()) {
		const size_t d = buffer_->tell() - pos;
		if (pos < buffer_->tell() && buffer_->begin + d <= buffer_->end) {
			buffer_->begin = buffer_->end - d;
			return *this;
		}
	}
	buffer_->seek(pos, SEEK_SET);
	return *this;
}

void Deserializer::seek_forward(size_t n)
{
	buffer_->seek(n, SEEK_CUR);
}

size_t Deserializer::read_raw(char *ptr, int64_t count)
{
	if (count <= avail()) {
		pop(ptr, count);
		return count;
	}
	size_t total = 0;
	do {
		const size_t n = std::min(count, avail());
		pop(ptr, n);
		ptr += n;
		count -= n;
		total += n;
		if (avail() == 0)
			buffer_->fetch();
	} while (count > 0 && avail() > 0);
	return total;
}

string Deserializer::peek(int64_t n) {
	if (avail() >= n)
		return string(buffer_->begin, n);
	if (avail() == 0)
		buffer_->fetch();
	if(avail() < n && !buffer_->eof())
		throw runtime_error("Invalid peek");
	return string(buffer_->begin, std::min(buffer_->begin + n, buffer_->end));
}

void Deserializer::pop(char *dst, int64_t n)
{
	assert(n <= avail());
	memcpy(dst, buffer_->begin, n);
	buffer_->begin += n;
}

bool Deserializer::seek_forward(char delimiter)
{
	int d = delimiter;
	do {
		const char *p = (const char*)memchr((void*)buffer_->begin, d, avail());
		if (p) {
			const size_t n = p - buffer_->begin;
			buffer_->begin += n + 1;
			return true;
		}
	} while (buffer_->fetch());
	return false;
}