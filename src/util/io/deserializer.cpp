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
#include <algorithm>
#include <string.h>
#include <iterator>
#include "deserializer.h"
#include "record_reader.h"

using namespace std;

Deserializer::Deserializer(StreamEntity* buffer):
	varint(false),
	buffer_(buffer),
	begin_(NULL),
	end_(NULL)
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
	begin_ = end_ = nullptr;
}

Deserializer& Deserializer::seek(int64_t pos)
{
	if (buffer_->seekable() && buffer_->tell()) {
		const size_t d = buffer_->tell() - pos;
		if (pos < buffer_->tell() && begin_ + d <= end_) {
			begin_ = end_ - d;
			return *this;
		}
	}
	buffer_->seek(pos, SEEK_SET);
	begin_ = end_ = nullptr;
	return *this;
}

void Deserializer::seek_forward(size_t n)
{
	buffer_->seek_forward(n);
	begin_ = end_ = nullptr;
}

size_t Deserializer::read_raw(char *ptr, size_t count)
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
			fetch();
	} while (count > 0 && avail() > 0);
	return total;
}

bool Deserializer::fetch()
{
	if (buffer_ == NULL)
		throw EndOfStream();
	pair<const char*, const char*> in = buffer_->read();
	begin_ = in.first;
	end_ = in.second;
	return end_ > begin_;
}

void Deserializer::pop(char *dst, size_t n)
{
	assert(n <= avail());
	memcpy(dst, begin_, n);
	begin_ += n;
}

bool Deserializer::seek_forward(char delimiter)
{
	int d = delimiter;
	do {
		const char *p = (const char*)memchr((void*)begin_, d, avail());
		if (p) {
			const size_t n = p - begin_;
			begin_ += n + 1;
			return true;
		}
	} while (fetch());
	return false;
}

DynamicRecordReader Deserializer::read_record() {
	return DynamicRecordReader(*this);
}