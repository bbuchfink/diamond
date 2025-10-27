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

#include <limits>
#include <stdexcept>
#include "compressed_buffer.h"
#ifdef WITH_ZSTD
#include <zstd.h>
#else
#include <zlib.h>
#endif

using std::runtime_error;

CompressedBuffer::CompressedBuffer() :
	buf_(BUF_SIZE)
{
	clear();
}

CompressedBuffer::~CompressedBuffer() {
	if (!stream_)
		return;
#ifdef WITH_ZSTD
	ZSTD_freeCStream(stream_);
#else
	deflateEnd(stream_);
	delete stream_;
#endif
}

void CompressedBuffer::write(const char* ptr, size_t n) {
#ifdef WITH_ZSTD
	ZSTD_inBuffer in_buf;
	in_buf.src = ptr;
	in_buf.size = n;
	in_buf.pos = 0;
	ZSTD_outBuffer out_buf;
	do {
		out_buf.dst = buf_.data();
		out_buf.size = buf_.size();
		out_buf.pos = size_;
		if (ZSTD_isError(ZSTD_compressStream(stream_, &out_buf, &in_buf)))
			throw runtime_error("ZSTD_compressStream");
		size_ = out_buf.pos;
		if (in_buf.pos < in_buf.size)
			buf_.resize(buf_.size() + BUF_SIZE);
	} while (in_buf.pos < in_buf.size);
#else
	const unsigned char* in = reinterpret_cast<const unsigned char*>(ptr);
	while (n > 0) {
		uInt chunk = static_cast<uInt>(std::min<int64_t>(n, std::numeric_limits<uInt>::max()));
		stream_->next_in = const_cast<Bytef*>(in);
		stream_->avail_in = chunk;
		while (stream_->avail_in > 0) {
			if (size_ == buf_.size())
				buf_.resize(buf_.size() + BUF_SIZE);
			uInt avail_out_before = static_cast<uInt>(buf_.size() - size_);
			stream_->next_out = reinterpret_cast<Bytef*>(buf_.data() + size_);
			stream_->avail_out = avail_out_before;
			int ret = deflate(stream_, Z_NO_FLUSH);
			if (ret == Z_STREAM_ERROR)
				throw std::runtime_error("deflate(Z_NO_FLUSH) failed");
			size_ += static_cast<size_t>(avail_out_before - stream_->avail_out);
			if (stream_->avail_in > 0 && stream_->avail_out == 0) {
				buf_.resize(buf_.size() + BUF_SIZE);
			}
		}
		in += chunk;
		n -= chunk;
	}
#endif
}

void CompressedBuffer::finish() {
#ifdef WITH_ZSTD
	ZSTD_outBuffer out_buf;
	size_t n;
	do {
		out_buf.dst = buf_.data();
		out_buf.size = buf_.size();
		out_buf.pos = size_;
		if (ZSTD_isError(n = ZSTD_endStream(stream_, &out_buf)))
			throw runtime_error("ZSTD_endStream");
		size_ = out_buf.pos;
		if (n > 0)
			buf_.resize(buf_.size() + BUF_SIZE);
	} while (n > 0);
	ZSTD_freeCStream(stream_);
#else
	int ret;
	stream_->next_in = Z_NULL;
	stream_->avail_in = 0;
	do {
		if (size_ == static_cast<size_t>(buf_.size()))
			buf_.resize(buf_.size() + BUF_SIZE);
		uInt avail_out_before = static_cast<uInt>(buf_.size() - size_);
		stream_->next_out = reinterpret_cast<Bytef*>(buf_.data() + size_);
		stream_->avail_out = avail_out_before;
		ret = deflate(stream_, Z_FINISH);
		if (ret != Z_OK && ret != Z_STREAM_END)
			throw runtime_error("deflate(Z_FINISH) failed");
		size_ += static_cast<size_t>(avail_out_before - stream_->avail_out);
		if (ret != Z_STREAM_END && stream_->avail_out == 0) {
			buf_.resize(buf_.size() + BUF_SIZE);
		}
	} while (ret != Z_STREAM_END);
	deflateEnd(stream_);
	delete stream_;
#endif
	stream_ = nullptr;
}

void CompressedBuffer::clear() {
#ifdef WITH_ZSTD
	stream_ = ZSTD_createCStream();
	if (!stream_)
		throw runtime_error("ZSTD_createCStream error");
#else
	stream_ = new z_stream;
	stream_->zalloc = Z_NULL;
	stream_->zfree = Z_NULL;
	stream_->opaque = Z_NULL;
	int ret = deflateInit(stream_, Z_DEFAULT_COMPRESSION);
	if (ret != Z_OK)
		throw runtime_error("deflateInit error");
#endif
	size_ = 0;
}