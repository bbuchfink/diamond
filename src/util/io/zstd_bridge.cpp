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
	if (!strm_) {
		strm_ = ZSTD_createDStream();
		if (!strm_)
			throw std::runtime_error("Error initializing zstd decompressor (ZSTD_createDStream)");
		in_.resize(chunk_size);
		out_.resize(chunk_size);
	}
	const size_t ret = ZSTD_initDStream(strm_);
	if (ZSTD_isError(ret))
		throw std::runtime_error(std::string("Error initializing zstd decompressor: ") + ZSTD_getErrorName(ret));
	out_pos_ = 0;
	out_end_ = 0;
	in_buf_ = { nullptr, 0, 0 };
	eos_ = false;
	pushback_ = EOF;
}

bool ZstdDecompressor::fill_out_buf(FILE* stream) {
	out_pos_ = 0;
	out_end_ = 0;
	ZSTD_outBuffer out_buf{ out_.data(), out_.size(), 0 };
	while (out_buf.pos == 0 && !eos_) {
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
	out_end_ = out_buf.pos;
	return out_end_ > 0;
}

size_t ZstdDecompressor::fread(void* buffer, size_t size, size_t count, FILE* stream) {
	if (size == 0 || count == 0)
		return 0;
	const size_t total = size * count;
	unsigned char* out = static_cast<unsigned char*>(buffer);
	size_t produced = 0;
	if (pushback_ != EOF) {
		out[produced++] = static_cast<unsigned char>(pushback_);
		pushback_ = EOF;
	}
	while (produced < total) {
		if (out_pos_ < out_end_) {
			const size_t n = std::min(total - produced, out_end_ - out_pos_);
			std::memcpy(out + produced, out_.data() + out_pos_, n);
			out_pos_ += n;
			produced += n;
		} else if (!fill_out_buf(stream)) {
			break;
		}
	}
	return produced / size;
}

int ZstdDecompressor::fgetc(FILE* stream) {
	if (pushback_ != EOF) {
		const int c = pushback_;
		pushback_ = EOF;
		return c;
	}
	if (out_pos_ < out_end_ || fill_out_buf(stream))
		return out_[out_pos_++];
	return EOF;
}

ssize_t ZstdDecompressor::getdelim(char** buf, size_t* buf_size, char delimiter, FILE* fp) {
	if (buf == nullptr || buf_size == nullptr)
		return -1;
	if (*buf == nullptr || *buf_size == 0) {
		*buf_size = 128;
		*buf = (char*)realloc(*buf, *buf_size);
		if (*buf == nullptr)
			throw std::runtime_error("Memory allocation failed in file I/O.");
	}
	size_t len = 0;
	if (pushback_ != EOF) {
		const int c = pushback_;
		pushback_ = EOF;
		if (len + 1 >= *buf_size) {
			const size_t new_size = *buf_size * 2;
			char* tmp = (char*)realloc(*buf, new_size);
			if (tmp == nullptr)
				throw std::runtime_error("Memory allocation failed in file I/O.");
			*buf = tmp;
			*buf_size = new_size;
		}
		(*buf)[len++] = (char)c;
		if ((char)c == delimiter) {
			(*buf)[len] = '\0';
			return (ssize_t)len;
		}
	}
	while (out_pos_ < out_end_ || fill_out_buf(fp)) {
		const unsigned char* start = out_.data() + out_pos_;
		const size_t avail = out_end_ - out_pos_;
		const unsigned char* found = static_cast<const unsigned char*>(std::memchr(start, (unsigned char)delimiter, avail));
		const size_t chunk = found ? static_cast<size_t>(found - start + 1) : avail;
		if (len + chunk + 1 > *buf_size) {
			size_t new_size = *buf_size;
			while (new_size < len + chunk + 1)
				new_size *= 2;
			char* tmp = (char*)realloc(*buf, new_size);
			if (tmp == nullptr)
				throw std::runtime_error("Memory allocation failed in file I/O.");
			*buf = tmp;
			*buf_size = new_size;
		}
		std::memcpy(*buf + len, start, chunk);
		len += chunk;
		out_pos_ += chunk;
		if (found != nullptr)
			break;
	}
	if (len == 0)
		return -1;
	(*buf)[len] = '\0';
	return (ssize_t)len;
}

size_t ZstdDecompressor::read_to_fasta_record_end(std::string& dst, FILE* stream) {
	size_t n = 0;
	if (pushback_ != EOF) {
		const int c = pushback_;
		if (c == '>' && !dst.empty() && (dst.back() == '\n' || dst.back() == '\r'))
			return 0;
		pushback_ = EOF;
		dst.push_back(static_cast<char>(c));
		++n;
	}
	while (out_pos_ < out_end_ || fill_out_buf(stream)) {
		const unsigned char* start = out_.data() + out_pos_;
		const size_t avail = out_end_ - out_pos_;
		const unsigned char* found = static_cast<const unsigned char*>(std::memchr(start, '>', avail));
		if (!found) {
			dst.append(reinterpret_cast<const char*>(start), avail);
			n += avail;
			out_pos_ = out_end_;
		} else {
			const size_t prefix = static_cast<size_t>(found - start);
			const char prev = prefix > 0 ? static_cast<char>(found[-1]) : (dst.empty() ? '\0' : dst.back());
			if (prev == '\n' || prev == '\r') {
				dst.append(reinterpret_cast<const char*>(start), prefix);
				n += prefix;
				out_pos_ += prefix + 1;
				pushback_ = '>';
				return n;
			}
			dst.append(reinterpret_cast<const char*>(start), prefix + 1);
			n += prefix + 1;
			out_pos_ += prefix + 1;
		}
	}
	return n;
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
