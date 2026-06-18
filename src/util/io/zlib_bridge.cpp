/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <string.h>
#include <errno.h>
#include <limits>
#include "compressor.h"

ZlibDecompressor::ZlibDecompressor() {
	reset();
}

void ZlibDecompressor::reset() {
	if (initialized_)
		inflateEnd(&strm_);
	strm_.zalloc = Z_NULL;
	strm_.zfree = Z_NULL;
	strm_.opaque = Z_NULL;
	strm_.avail_in = 0;
	strm_.next_in = Z_NULL;
	if (inflateInit2(&strm_, 15 + 32) != Z_OK)
		throw std::runtime_error("Error initializing zlib decompressor (inflateInit2)");
	in_.resize(chunk_size);
	eos_ = false;
	initialized_ = true;
}

size_t ZlibDecompressor::fread(void* buffer, size_t size, size_t count, FILE* stream) {
	if (size == 0 || count == 0)
		return 0;
	const size_t total = size * count;
	unsigned char* out = static_cast<unsigned char*>(buffer);
	size_t produced = 0;
	if (pushback_ != EOF) {
		out[produced++] = static_cast<unsigned char>(pushback_);
		pushback_ = EOF;
	}
	while (produced < total && !eos_) {
		if (strm_.avail_in == 0) {
			const size_t rd = std::fread(in_.data(), 1, in_.size(), stream);
			if (rd == 0) {
				if (std::ferror(stream))
					throw std::runtime_error(std::string("Error reading compressed file: ") + strerror(errno));
				eos_ = true;
				break;
			}
			strm_.next_in = reinterpret_cast<Bytef*>(in_.data());
			strm_.avail_in = static_cast<uInt>(rd);
		}

		const uInt out_chunk = static_cast<uInt>(std::min(total - produced, (size_t)std::numeric_limits<uInt>::max()));
		strm_.next_out = reinterpret_cast<Bytef*>(out + produced);
		strm_.avail_out = out_chunk;

		const int ret = inflate(&strm_, Z_NO_FLUSH);
		produced += out_chunk - strm_.avail_out;

		if (ret == Z_STREAM_END) {
#if ZLIB_VERNUM >= 0x1270
			if (inflateReset2(&strm_, 15 + 32) != Z_OK)
				throw std::runtime_error("Error resetting compressed stream (inflateReset2)");
#else
			if (inflateReset(&strm_) != Z_OK)
				throw std::runtime_error("Error resetting compressed stream (inflateReset)");
#endif
			continue;
		}
		if (ret == Z_BUF_ERROR) {
			continue;
		}
		if (ret != Z_OK)
			throw std::runtime_error(std::string("Error during zlib decompression. The file may be corrupted: ") + (strm_.msg ? strm_.msg : "unknown error"));
	}
	return produced / size;
}

int ZlibDecompressor::fgetc(FILE* stream) {
	unsigned char c;
	return fread(&c, 1, 1, stream) == 1 ? c : EOF;
}

ssize_t ZlibDecompressor::getdelim(char** buf, size_t* buf_size, char delimiter, FILE* fp) {
	return getdelim_generic(buf, buf_size, delimiter, [this, fp] { return this->fgetc(fp); });
}

int ZlibDecompressor::ungetc(int c, FILE* stream) {
	if (c == EOF || pushback_ != EOF)
		return EOF;
	pushback_ = static_cast<unsigned char>(c);
	return c;
}

ZlibDecompressor::~ZlibDecompressor() {
	if (initialized_)
		inflateEnd(&strm_);
}

ZlibCompressor::ZlibCompressor() {
	strm_.zalloc = Z_NULL;
	strm_.zfree = Z_NULL;
	strm_.opaque = Z_NULL;
	strm_.avail_in = 0;
	strm_.next_in = Z_NULL;
	if (deflateInit2(&strm_, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK)
		throw std::runtime_error("Error initializing zlib compressor (deflateInit2)");
	out_.resize(chunk_size);
	initialized_ = true;
}

void ZlibCompressor::deflate_loop(FILE* stream, int flush) {
	int ret;
	do {
		strm_.next_out = reinterpret_cast<Bytef*>(out_.data());
		strm_.avail_out = static_cast<uInt>(out_.size());
		ret = deflate(&strm_, flush);
		if (ret == Z_STREAM_ERROR)
			throw std::runtime_error("Error during zlib compression (deflate)");
		const size_t produced = out_.size() - strm_.avail_out;
		if (produced != 0 && std::fwrite(out_.data(), 1, produced, stream) != produced)
			throw std::runtime_error(std::string("Error writing compressed file: ") + strerror(errno));
	} while (flush == Z_FINISH ? ret != Z_STREAM_END : strm_.avail_in != 0 || strm_.avail_out == 0);
}

size_t ZlibCompressor::fwrite(const void* buffer, size_t size, size_t count, FILE* stream) {
	if (size == 0 || count == 0)
		return 0;
	if (closed_)
		throw std::runtime_error("Cannot write to closed zlib compressor.");
	if (stream_ == nullptr)
		stream_ = stream;
	else if (stream_ != stream)
		throw std::runtime_error("Cannot write one zlib stream to multiple FILE handles.");
	if (count > std::numeric_limits<size_t>::max() / size)
		throw std::runtime_error("Compressed write size overflow.");

	size_t remaining = size * count;
	const unsigned char* in = static_cast<const unsigned char*>(buffer);
	while (remaining > 0) {
		const uInt chunk = static_cast<uInt>(std::min(remaining, (size_t)std::numeric_limits<uInt>::max()));
		strm_.next_in = (Bytef*)(in);
		strm_.avail_in = chunk;
		deflate_loop(stream, Z_NO_FLUSH);
		in += chunk;
		remaining -= chunk;
	}
	return count;
}

void ZlibCompressor::close(FILE* stream) {
	if (closed_)
		return;
	if (stream_ == nullptr)
		stream_ = stream;
	else if (stream_ != stream)
		throw std::runtime_error("Cannot close zlib compressor with a different FILE handle.");

	try {
		strm_.next_in = Z_NULL;
		strm_.avail_in = 0;
		deflate_loop(stream_, Z_FINISH);
	}
	catch (...) {
		deflateEnd(&strm_);
		initialized_ = false;
		closed_ = true;
		throw;
	}
	deflateEnd(&strm_);
	initialized_ = false;
	closed_ = true;
}

ZlibCompressor::~ZlibCompressor() {
	if (!initialized_)
		return;
	if (!closed_ && stream_) {
		try {
			strm_.next_in = Z_NULL;
			strm_.avail_in = 0;
			deflate_loop(stream_, Z_FINISH);
		}
		catch (...) {
		}
	}
	deflateEnd(&strm_);
}
