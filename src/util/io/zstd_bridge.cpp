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

#ifdef WITH_ZSTD

ZstdDecompressor::ZstdDecompressor() {
	reset();
}

void ZstdDecompressor::reset() {
	if (strm_)
		ZSTD_freeDStream(strm_);
	strm_ = ZSTD_createDStream();
	if (!strm_)
		throw std::runtime_error("Error initializing zstd decompressor (ZSTD_createDStream)");
	const size_t ret = ZSTD_initDStream(strm_);
	if (ZSTD_isError(ret))
		throw std::runtime_error(std::string("Error initializing zstd decompressor: ") + ZSTD_getErrorName(ret));
	in_.resize(chunk_size);
	eos_ = false;
}

size_t ZstdDecompressor::fread(void* buffer, size_t size, size_t count, FILE* stream) {
	if (size == 0 || count == 0)
		return 0;
	const size_t total = size * count;
	unsigned char* out = static_cast<unsigned char*>(buffer);
	ZSTD_outBuffer out_buf{ out, total, 0 };
	if (pushback_ != EOF) {
		out[out_buf.pos++] = static_cast<unsigned char>(pushback_);
		pushback_ = EOF;
	}
	while (out_buf.pos < total && !eos_) {
		if (in_buf_.pos == in_buf_.size) {
			const size_t rd = std::fread(in_.data(), 1, in_.size(), stream);
			if (rd == 0) {
				if (std::ferror(stream))
					throw std::runtime_error(std::string("Error reading compressed file: ") + strerror(errno));
				eos_ = true;
				break;
			}
			in_buf_.src = in_.data();
			in_buf_.size = rd;
			in_buf_.pos = 0;
		}
		const size_t ret = ZSTD_decompressStream(strm_, &out_buf, &in_buf_);
		if (ZSTD_isError(ret))
			throw std::runtime_error(std::string("Error during zstd decompression. The file may be corrupted: ") + ZSTD_getErrorName(ret));
	}
	return out_buf.pos / size;
}

int ZstdDecompressor::fgetc(FILE* stream) {
	unsigned char c;
	return fread(&c, 1, 1, stream) == 1 ? c : EOF;
}

ssize_t ZstdDecompressor::getdelim(char** buf, size_t* buf_size, char delimiter, FILE* fp) {
	return getdelim_generic(buf, buf_size, delimiter, [this, fp] { return this->fgetc(fp); });
}

int ZstdDecompressor::ungetc(int c, FILE* stream) {
	if (c == EOF || pushback_ != EOF)
		return EOF;
	pushback_ = static_cast<unsigned char>(c);
	return c;
}

ZstdDecompressor::~ZstdDecompressor() {
	if (strm_)
		ZSTD_freeDStream(strm_);
}

ZstdCompressor::ZstdCompressor() {
	strm_ = ZSTD_createCStream();
	if (!strm_)
		throw std::runtime_error("Error initializing zstd compressor (ZSTD_createCStream)");
	const size_t ret = ZSTD_initCStream(strm_, 0);
	if (ZSTD_isError(ret))
		throw std::runtime_error(std::string("Error initializing zstd compressor: ") + ZSTD_getErrorName(ret));
	out_.resize(chunk_size);
}

void ZstdCompressor::write_out(FILE* stream, size_t count) {
	if (count != 0 && std::fwrite(out_.data(), 1, count, stream) != count)
		throw std::runtime_error(std::string("Error writing compressed file: ") + strerror(errno));
}

size_t ZstdCompressor::fwrite(const void* buffer, size_t size, size_t count, FILE* stream) {
	if (size == 0 || count == 0)
		return 0;
	if (closed_)
		throw std::runtime_error("Cannot write to closed zstd compressor.");
	if (stream_ == nullptr)
		stream_ = stream;
	else if (stream_ != stream)
		throw std::runtime_error("Cannot write one zstd stream to multiple FILE handles.");
	if (count > std::numeric_limits<size_t>::max() / size)
		throw std::runtime_error("Compressed write size overflow.");

	ZSTD_inBuffer in_buf{ buffer, size * count, 0 };
	while (in_buf.pos < in_buf.size) {
		ZSTD_outBuffer out_buf{ out_.data(), out_.size(), 0 };
		const size_t ret = ZSTD_compressStream(strm_, &out_buf, &in_buf);
		if (ZSTD_isError(ret))
			throw std::runtime_error(std::string("Error during zstd compression: ") + ZSTD_getErrorName(ret));
		write_out(stream, out_buf.pos);
	}
	return count;
}

void ZstdCompressor::close(FILE* stream) {
	if (closed_)
		return;
	if (stream_ == nullptr)
		stream_ = stream;
	else if (stream_ != stream)
		throw std::runtime_error("Cannot close zstd compressor with a different FILE handle.");

	try {
		size_t remaining;
		do {
			ZSTD_outBuffer out_buf{ out_.data(), out_.size(), 0 };
			remaining = ZSTD_endStream(strm_, &out_buf);
			if (ZSTD_isError(remaining))
				throw std::runtime_error(std::string("Error finishing zstd compression: ") + ZSTD_getErrorName(remaining));
			write_out(stream_, out_buf.pos);
		} while (remaining != 0);
	}
	catch (...) {
		ZSTD_freeCStream(strm_);
		strm_ = nullptr;
		closed_ = true;
		throw;
	}
	ZSTD_freeCStream(strm_);
	strm_ = nullptr;
	closed_ = true;
}

ZstdCompressor::~ZstdCompressor() {
	if (!strm_)
		return;
	if (!closed_ && stream_) {
		try {
			size_t remaining;
			do {
				ZSTD_outBuffer out_buf{ out_.data(), out_.size(), 0 };
				remaining = ZSTD_endStream(strm_, &out_buf);
				if (ZSTD_isError(remaining))
					break;
				write_out(stream_, out_buf.pos);
			} while (remaining != 0);
		}
		catch (...) {
		}
	}
	ZSTD_freeCStream(strm_);
}

#endif
