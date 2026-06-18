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
#include "decompressor.h"

struct CompressorX {
	virtual size_t fwrite(const void* buffer, size_t size, size_t count, FILE* stream) = 0;	
	virtual void close(FILE* stream) {}
	virtual CompressionLib lib() const = 0;
	virtual ~CompressorX() = default;
};

struct PassThroughCompressor : CompressorX {
	size_t fwrite(const void* buffer, size_t size, size_t count, FILE* stream) override {
		return std::fwrite(buffer, size, count, stream);
	}	
	virtual CompressionLib lib() const override {
		return CompressionLib::NONE;
	}
};

struct ZlibCompressor : CompressorX {
	ZlibCompressor();
	size_t fwrite(const void* buffer, size_t size, size_t count, FILE* stream) override;
	void close(FILE* stream) override;
	virtual CompressionLib lib() const override {
		return CompressionLib::ZLIB;
	}
	~ZlibCompressor();
private:
	void deflate_loop(FILE* stream, int flush);
	static const size_t chunk_size = 1llu << 20;
	z_stream strm_;
	std::vector<unsigned char> out_;
	FILE* stream_ = nullptr;
	bool initialized_ = false;
	bool closed_ = false;
};

#ifdef WITH_ZSTD
struct ZstdCompressor : CompressorX {
	ZstdCompressor();
	size_t fwrite(const void* buffer, size_t size, size_t count, FILE* stream) override;
	void close(FILE* stream) override;
	virtual CompressionLib lib() const override {
		return CompressionLib::ZSTD;
	}
	~ZstdCompressor();
private:
	void write_out(FILE* stream, size_t count);
	static const size_t chunk_size = 1llu << 20;
	ZSTD_CStream* strm_ = nullptr;
	std::vector<unsigned char> out_;
	FILE* stream_ = nullptr;
	bool closed_ = false;
};
#endif
