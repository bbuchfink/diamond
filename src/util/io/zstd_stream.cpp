/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#include <vector>
#include <string.h>
#include <errno.h>
#include "zstd_stream.h"
#include "../system.h"

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
	size_t zr = ZSTD_initDStream(stream);
	if (ZSTD_isError(zr)) throw runtime_error(ZSTD_getErrorName(zr));
}

size_t ZstdSource::read(char* ptr, size_t count)
{
	ZSTD_outBuffer out_buf{ ptr, count, 0 };
	auto* buf = static_cast<InputStreamBuffer*>(prev_);

	while (out_buf.pos < out_buf.size) {
		if (in_buf.pos == in_buf.size) {
			if (buf->begin == buf->end)
				buf->fetch();
			in_buf.src = buf->begin;
			in_buf.size = static_cast<size_t>(buf->end - buf->begin);
			in_buf.pos = 0;

			if (in_buf.size == 0) {
				eos_ = true;
				break;
			}
		}

		size_t const ret = ZSTD_decompressStream(stream, &out_buf, &in_buf);
		if (ZSTD_isError(ret))
			throw runtime_error(string("ZSTD_decompressStream: ")
				+ ZSTD_getErrorName(ret));
		size_t const consumed = in_buf.pos;
		buf->begin += consumed;
		size_t const remaining = in_buf.size - consumed;
		if (remaining) {
			in_buf.src = buf->begin;
			in_buf.size = remaining;
			in_buf.pos = 0;
		}
		else {
			in_buf.src = nullptr;
			in_buf.size = 0;
			in_buf.pos = 0;
			//if (buf->begin == buf->end && ret > 0) {
				//throw runtime_error("truncated zstd stream (need more input)");
		}
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

size_t zstd_decompress(FILE* src, void* dst, size_t dstCapacity) {
	ZSTD_DCtx* dctx = ZSTD_createDCtx();
	if (!dctx)
		throw std::runtime_error("ZSTD_createDCtx");
	const size_t inChunkSize = ZSTD_DStreamInSize();
	std::vector<unsigned char> in(inChunkSize);
	unsigned char* outBase = static_cast<unsigned char*>(dst);
	size_t totalOut = 0;
	size_t lastRet = 1;
	for (;;) {
		const size_t read = std::fread(in.data(), 1, in.size(), src);
		if (std::ferror(src)) {
			ZSTD_freeDCtx(dctx);
			throw std::runtime_error(string("Error reading file: ") + strerror(errno));
		}
		ZSTD_inBuffer input{ in.data(), read, 0 };
		while (input.pos < input.size) {
			ZSTD_outBuffer output{ outBase + totalOut, dstCapacity - totalOut, 0 };
			size_t const ret = ZSTD_decompressStream(dctx, &output, &input);
			if (ZSTD_isError(ret)) {
				ZSTD_freeDCtx(dctx);
				throw std::runtime_error(ZSTD_getErrorName(ret));
			}
			totalOut += output.pos;
			lastRet = ret;
			if (totalOut > dstCapacity) {
				ZSTD_freeDCtx(dctx);
				throw std::runtime_error("Failed decompressing zstd stream: output buffer too small");
			}
		}
		if (read == 0) {
			if (lastRet != 0) {
				ZSTD_freeDCtx(dctx);
				throw std::runtime_error("Failed decompressing zstd stream");
			}
			break;
		}
	}
	ZSTD_freeDCtx(dctx);
	return totalOut;
}