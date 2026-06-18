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

#pragma once
#include <vector>
#include <stdexcept>
#include <zlib.h>
#ifdef WITH_ZSTD
#include <zstd.h>
#endif

#ifdef _MSC_VER
using ssize_t = int64_t;
#endif

enum class CompressionLib { NONE, ZLIB, ZSTD };

struct Decompressor {
	virtual size_t fread(void* buffer, size_t size, size_t count, FILE* stream) = 0;
	virtual int fgetc(FILE* stream) = 0;
	virtual ssize_t getdelim(char** buf, size_t* buf_size, char delimiter, FILE* fp) = 0;
	ssize_t getline(char** buf, size_t* buf_size, FILE* fp) {
		return getdelim(buf, buf_size, '\n', fp);
	}
	virtual CompressionLib lib() const = 0;
	virtual int ungetc(int c, FILE* stream) = 0;
	virtual void reset() = 0;
};

template<typename ReadByte>
inline ssize_t getdelim_generic(char** buf, size_t* buf_size, char delimiter, ReadByte read_byte) {
	if (buf == nullptr || buf_size == nullptr)
		return -1;
	if (*buf == nullptr || *buf_size == 0) {
		*buf_size = 128;
		*buf = (char*)realloc(*buf, *buf_size);
		if (*buf == nullptr)
			throw std::runtime_error("Memory allocation failed in file I/O.");
	}
	ssize_t len = 0;
	int c;
	while ((c = read_byte()) != EOF) {
		if ((size_t)len + 1 >= *buf_size) {
			const size_t new_size = *buf_size * 2;
			char* tmp = (char*)realloc(*buf, new_size);
			if (tmp == nullptr)
				throw std::runtime_error("Memory allocation failed in file I/O.");
			*buf = tmp;
			*buf_size = new_size;
		}
		(*buf)[len++] = (char)c;
		if ((char)c == delimiter)
			break;
	}
	if (len == 0)
		return -1;
	(*buf)[len] = '\0';
	return len;
}

struct PassThrough : Decompressor {
	size_t fread(void* buffer, size_t size, size_t count, FILE* stream) override {
		return std::fread(buffer, size, count, stream);
	}
	int fgetc(FILE* stream) override {
#if defined(_MSC_VER) || defined(__APPLE__)
		return std::fgetc(stream);
#else
		return fgetc_unlocked(stream);
#endif
	}
	ssize_t getdelim(char** buf, size_t* buf_size, char delimiter, FILE* fp) override {
#ifdef _MSC_VER
		return getdelim_generic(buf, buf_size, delimiter, [fp] { return std::fgetc(fp); });
#else
		return ::getdelim(buf, buf_size, delimiter, fp);
#endif
	}
	virtual CompressionLib lib() const override {
		return CompressionLib::NONE;
	}
	virtual int ungetc(int c, FILE* stream) override {
		return std::ungetc(c, stream);
	}
	virtual void reset() override {}
};

struct ZlibDecompressor : Decompressor {
	ZlibDecompressor();
	size_t fread(void* buffer, size_t size, size_t count, FILE* stream) override;
	int fgetc(FILE* stream) override;
	ssize_t getdelim(char** buf, size_t* buf_size, char delimiter, FILE* fp) override;
	virtual CompressionLib lib() const override {
		return CompressionLib::ZLIB;
	}
	int ungetc(int c, FILE* stream) override;
	virtual void reset() override;
	~ZlibDecompressor();
private:
	static const size_t chunk_size = 1llu << 20;
	z_stream strm_;
	std::vector<unsigned char> in_;
	bool eos_ = false;
	bool initialized_ = false;
	int pushback_ = EOF;
};

#ifdef WITH_ZSTD
struct ZstdDecompressor : Decompressor {
	ZstdDecompressor();
	size_t fread(void* buffer, size_t size, size_t count, FILE* stream) override;
	int fgetc(FILE* stream) override;
	ssize_t getdelim(char** buf, size_t* buf_size, char delimiter, FILE* fp) override;
	virtual CompressionLib lib() const override {
		return CompressionLib::ZSTD;
	}
	int ungetc(int c, FILE* stream) override;
	virtual void reset() override;
	~ZstdDecompressor();
private:
	static const size_t chunk_size = 1llu << 20;
	ZSTD_DStream* strm_ = nullptr;
	std::vector<unsigned char> in_;
	ZSTD_inBuffer in_buf_{ nullptr, 0, 0 };
	bool eos_ = false;
	int pushback_ = EOF;
};
#endif