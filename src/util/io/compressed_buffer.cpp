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
	ZSTD_inBuffer in_buf{ ptr, n, 0 };
	while (in_buf.pos < in_buf.size) {
		if (size_ == buf_.size())
			buf_.resize(buf_.size() + BUF_SIZE);
		ZSTD_outBuffer out_buf;
		out_buf.dst = buf_.data() + size_;
		out_buf.size = buf_.size() - size_;
		out_buf.pos = 0;
		size_t zr = ZSTD_compressStream(static_cast<ZSTD_CStream*>(stream_), &out_buf, &in_buf);
		if (ZSTD_isError(zr))
			throw runtime_error(ZSTD_getErrorName(zr));
		size_ += out_buf.pos;
		if (out_buf.pos == 0 && in_buf.pos < in_buf.size)
			buf_.resize(buf_.size() + BUF_SIZE);
	}
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
	for (;;) {
		if (size_ == buf_.size())
			buf_.resize(buf_.size() + BUF_SIZE);
		ZSTD_outBuffer out_buf;
		out_buf.dst = buf_.data() + size_;
		out_buf.size = buf_.size() - size_;
		out_buf.pos = 0;
		size_t remaining = ZSTD_endStream(static_cast<ZSTD_CStream*>(stream_), &out_buf);
		if (ZSTD_isError(remaining))
			throw runtime_error(ZSTD_getErrorName(remaining));
		size_ += out_buf.pos;
		if (remaining == 0) break;
	}
	ZSTD_freeCStream(static_cast<ZSTD_CStream*>(stream_));
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
	size_t zr = ZSTD_initCStream(stream_, 0);
	if (ZSTD_isError(zr)) throw runtime_error(ZSTD_getErrorName(zr));
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