/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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

FILE* Deserializer::file()
{
	return buffer_->file();
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