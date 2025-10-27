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

#include "zstd_stream.h"

using std::pair;
using std::runtime_error;
using std::string;

ZstdSink::ZstdSink(StreamEntity* prev) :
	StreamEntity(prev),
	stream(ZSTD_createCStream())
{
	if (!stream)
		throw runtime_error("ZSTD_createCStream error");
}

void ZstdSink::write(const char* ptr, size_t count)
{
	ZSTD_inBuffer in_buf;
	in_buf.src = ptr;
	in_buf.size = count;
	in_buf.pos = 0;
	ZSTD_outBuffer out_buf;
	do {
		pair<char*, char*> out = prev_->write_buffer();
		out_buf.dst = out.first;
		out_buf.size = out.second - out.first;
		out_buf.pos = 0;
		if(ZSTD_isError(ZSTD_compressStream(stream, &out_buf, &in_buf)))
			throw runtime_error("ZSTD_compressStream");
		prev_->flush(out_buf.pos);
	} while (in_buf.pos < in_buf.size);
}

void ZstdSink::close()
{
	if (!stream)
		return;
	ZSTD_outBuffer out_buf;
	size_t n;
	do {
		pair<char*, char*> out = prev_->write_buffer();
		out_buf.dst = out.first;
		out_buf.size = out.second - out.first;
		out_buf.pos = 0;
		if (ZSTD_isError(n = ZSTD_endStream(stream, &out_buf)))
			throw runtime_error("ZSTD_endStream");
		prev_->flush(out_buf.pos);
	} while (n > 0);
	ZSTD_freeCStream(stream);
	stream = nullptr;
	prev_->close();
}

ZstdSource::ZstdSource(InputStreamBuffer* prev):
StreamEntity(prev)
{
	init();
}


void ZstdSource::init()
{
	eos_ = false;
	stream = ZSTD_createDStream();
	in_buf.pos = 0;
	in_buf.size = 0;
	in_buf.src = nullptr;
}

size_t ZstdSource::read(char* ptr, size_t count)
{
	ZSTD_outBuffer out_buf;
	out_buf.dst = ptr;
	out_buf.pos = 0;
	out_buf.size = count;
	InputStreamBuffer* buf = static_cast<InputStreamBuffer*>(prev_);
	while (out_buf.pos < out_buf.size && !eos_) {
		if (in_buf.pos >= in_buf.size) {
			if (buf->begin == buf->end)
				buf->fetch();
			in_buf.pos = 0;
			in_buf.size = buf->end - buf->begin;
			in_buf.src = buf->begin;
			buf->begin += in_buf.size;
			if (in_buf.size == 0) {
				eos_ = true;
				break;
			}
		}
		const size_t code = ZSTD_decompressStream(stream, &out_buf, &in_buf);
		if (ZSTD_isError(code))
			throw runtime_error(string("ZSTD_decompressStream: ") + ZSTD_getErrorName(code));;
	}
	return out_buf.pos;
}

bool ZstdSource::eof() {
	return eos_;
}

void ZstdSource::close()
{
	if (!stream)
		return;
	ZSTD_freeDStream(stream);
	stream = nullptr;
	prev_->close();
}

void ZstdSource::rewind()
{
	prev_->rewind();
	init();
}