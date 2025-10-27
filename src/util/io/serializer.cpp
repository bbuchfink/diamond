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

#include <string.h>
#include "serializer.h"

using std::string;
using std::pair;

Serializer::Serializer(StreamEntity *buffer) :
	buffer_(buffer)
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

void Serializer::consume(const char *ptr, size_t n) {
	write_raw(ptr, n);
}

void Serializer::finalize() {
	close();
}