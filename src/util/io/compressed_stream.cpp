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

#include <stdexcept>
#include <cstdio>
#include <cerrno>
#include <errno.h>
#include <vector>
#include <limits>
#include <zlib.h>
#include <string.h>
#include "compressed_stream.h"
#include "../system.h"

using std::pair;
using std::string;

void ZlibSource::init()
{
	eos_ = false;
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	int ret = inflateInit2(&strm, 15 + 32);
	if (ret != Z_OK)
		throw std::runtime_error("Error opening compressed file (inflateInit): " + file_name());
}

ZlibSource::ZlibSource(InputStreamBuffer *prev):
	StreamEntity(prev)
{
	init();
}

size_t ZlibSource::read(char *ptr, size_t count)
{
	InputStreamBuffer* buf = static_cast<InputStreamBuffer*>(prev_);
	strm.avail_out = (uInt)count;
	strm.next_out = (Bytef*)ptr;
	while (strm.avail_out > 0 && !eos_) {
		if (strm.avail_in == 0) {
			if (buf->end == buf->begin)
				buf->fetch();

			strm.avail_in = (uInt)(buf->end - buf->begin);
			if (strm.avail_in == 0) {
				eos_ = true;
				break;
			}
			strm.next_in = (Bytef*)buf->begin;
			buf->begin += strm.avail_in;
		}

		int ret = inflate(&strm, Z_NO_FLUSH);
		if (ret == Z_STREAM_END) {
			int ret = inflateInit2(&strm, 15 + 32);
			if (ret != Z_OK)
				throw std::runtime_error("Error initializing compressed stream (inflateInit): " + file_name());
		}
		else if (ret != Z_OK)
			throw std::runtime_error("Error reading gzip-compressed input file. The file may be corrupted.");
	}
	return count - strm.avail_out;
}

bool ZlibSource::eof() {
	return eos_;
}

void ZlibSource::close()
{
	inflateEnd(&strm);
	prev_->close();
}

void ZlibSource::rewind()
{
	prev_->rewind();
	inflateEnd(&strm);
	init();
}

ZlibSink::ZlibSink(StreamEntity *prev):
	StreamEntity(prev)
{
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK)
		throw std::runtime_error("deflateInit error");
}

void ZlibSink::deflate_loop(const char *ptr, size_t count, int flush)
{
	strm.avail_in = (uInt)count;
	strm.next_in = (Bytef*)ptr;
	do {
		pair<char*, char*> out = prev_->write_buffer();
		const size_t chunk_size = out.second - out.first;
		strm.avail_out = (uInt)chunk_size;
		strm.next_out = (Bytef*)out.first;
		if (deflate(&strm, flush) == Z_STREAM_ERROR)
			throw std::runtime_error("deflate error");
		prev_->flush(chunk_size - strm.avail_out);
	} while (strm.avail_out == 0);
}

void ZlibSink::write(const char * ptr, size_t count)
{
	deflate_loop(ptr, count, Z_NO_FLUSH);
}

void ZlibSink::close()
{
	deflate_loop(0, 0, Z_FINISH);
	deflateEnd(&strm);
	prev_->close();
}

size_t zlib_decompress(FILE* src, void* dst, size_t dstCapacity) {
    z_stream zs{};
    int ret = inflateInit2(&zs, 15 + 32);	
    if (ret != Z_OK)
		throw std::runtime_error("inflateInit2");
    const size_t inChunkSize = 64 * 1024;
    std::vector<unsigned char> in(inChunkSize);
    unsigned char* outBase = static_cast<unsigned char*>(dst);
    size_t totalOut = 0;
	bool endedLastStream = false;
	for (;;) {
		size_t rd = std::fread(in.data(), 1, in.size(), src);
		if (std::ferror(src))
			throw std::runtime_error(string("Error reading file: ") + strerror(errno));
		zs.next_in = in.data();
		zs.avail_in = static_cast<uInt>(rd);
		while (zs.avail_in > 0) {
			uInt outChunk = (uInt)std::min(dstCapacity - totalOut, (size_t)std::numeric_limits<uInt>::max());
			if (outChunk == 0) {
				unsigned char scratch;
				zs.next_out = &scratch;
				zs.avail_out = 1;
				int z = inflate(&zs, Z_NO_FLUSH);
				if (z == Z_STREAM_END) {
					endedLastStream = true;
#if ZLIB_VERNUM >= 0x1270
					if (inflateReset2(&zs, 15 + 32) != Z_OK) { inflateEnd(&zs); throw std::runtime_error("inflateReset2"); }
#else
					if (inflateReset(&zs) != Z_OK) { inflateEnd(&zs); hard_fail("inflateReset"); }
#endif
					continue;
				}
				if (z == Z_OK && zs.avail_out == 0) {
					inflateEnd(&zs);
					throw std::runtime_error("zlib_decompress: output buffer too small");
				}
				if (z == Z_BUF_ERROR && zs.avail_in == 0) {
					break;
				}
				inflateEnd(&zs);
				throw std::runtime_error(std::string("Error during zlib decompression: ") + (zs.msg ? zs.msg : "unknown error"));
			}
			zs.next_out = outBase + totalOut;
			zs.avail_out = outChunk;
			int z = inflate(&zs, Z_NO_FLUSH);
			totalOut += (outChunk - zs.avail_out);
			if (z == Z_STREAM_END) {
				endedLastStream = true;
#if ZLIB_VERNUM >= 0x1270
				if (inflateReset2(&zs, 15 + 32) != Z_OK) { inflateEnd(&zs); throw std::runtime_error("inflateReset2"); }
#else
				if (inflateReset(&zs) != Z_OK) { inflateEnd(&zs); throw std::runtime_error("inflateReset"); }
#endif
				continue;
			}
			if (z == Z_OK) {
				endedLastStream = false;
				if (zs.avail_out == 0) continue;
				break;
			}
			if (z == Z_BUF_ERROR) {
				if (zs.avail_in == 0) break;
				if (zs.avail_out == 0) {
					inflateEnd(&zs);
					throw std::runtime_error("zlib_decompress: output buffer too small");
				}
			}
			if (z != Z_OK && z != Z_BUF_ERROR) {
				inflateEnd(&zs);
				throw std::runtime_error("Error during zlib decompression: " + std::string(zs.msg ? zs.msg : "unknown error"));
			}
		}
		if (rd == 0) {
			if (!endedLastStream) {
				inflateEnd(&zs);
				throw std::runtime_error("Unexpected end of zlib stream");
			}
			break;
		}
	}
	inflateEnd(&zs);
	return totalOut;
}