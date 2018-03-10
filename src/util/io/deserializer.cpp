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

#include <assert.h>
#include <algorithm>
#include <string.h>
#include "deserializer.h"

template<typename _t>
void append(_t &container, const char *ptr, size_t n)
{
	container.append(ptr, n);
}

template<>
void append<vector<char> >(vector<char> &container, const char *ptr, size_t n)
{
	container.insert(container.end(), ptr, ptr + n);
}

Deserializer::Deserializer(StreamEntity* buffer):
	buffer_(buffer),
	begin_(NULL),
	end_(NULL),
	varint(false)
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
	begin_ = end_ = NULL;
}

Deserializer& Deserializer::seek(size_t pos)
{
	buffer_->seek(pos);
	begin_ = end_ = NULL;
	return *this;
}

void Deserializer::seek_forward(size_t n)
{
	buffer_->seek_forward(n);
	begin_ = end_ = NULL;
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

template<typename _t>
bool Deserializer::read_to(_t &container, char delimiter)
{
	int d = delimiter;
	container.clear();
	do {
		const char *p = (const char*)memchr((void*)begin_, d, avail());
		if (p == 0)
			append(container, begin_, avail());
		else {
			const size_t n = p - begin_;
			append(container, begin_, n);
			begin_ += n + 1;
			return true;
		}
	} while (fetch());
	return !container.empty();
}

bool Deserializer::read_until(string &dst, char delimiter)
{
	return read_to(dst, delimiter);
}

bool Deserializer::read_until(vector<char> &dst, char delimiter)
{
	return read_to(dst, delimiter);
}