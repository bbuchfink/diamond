/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
						Benjamin Buchfink

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include <iostream>
#ifdef __APPLE__
#include <TargetConditionals.h>
#elif not defined(_MSC_VER)
#include <fcntl.h>
#endif
#include "output_file.h"
#include "file_sink.h"
#include "output_stream_buffer.h"
#include "compressed_stream.h"
#ifdef WITH_ZSTD
#include "zstd_stream.h"
#endif
#include "../system.h"

using std::string;

static StreamEntity* make_compressor(const Compressor c, StreamEntity* buffer) {
	switch (c) {
	case Compressor::ZLIB:
		return new ZlibSink(buffer);
	case Compressor::ZSTD:
#ifdef WITH_ZSTD
		return new ZstdSink(buffer);
#else
		throw std::runtime_error("Executable was not compiled with ZStd support.");
#endif
	default:
		throw std::runtime_error("");
	}
}

OutputFile::OutputFile(const string &file_name, Compressor compressor, const char *mode) :
	Serializer(new OutputStreamBuffer(new FileSink(file_name, mode))),
	file_name_(file_name)
{
	if (compressor != Compressor::NONE) {
		buffer_ = new OutputStreamBuffer(make_compressor(compressor, buffer_));
		reset_buffer();
	}
}

OutputFile::OutputFile(const TempFileData& d, Compressor compressor, const char *mode):
#ifdef _MSC_VER
	OutputFile(d.name, compressor, mode)
#else
	Serializer(new OutputStreamBuffer(new FileSink(d.name, d.fd, mode))),
	file_name_(d.name)
#endif	
{
}

void OutputFile::remove()
{
	if (::remove(file_name_.c_str()) != 0)
		std::cerr << "Warning: Failed to delete file " << file_name_ << std::endl;
}

void OutputFile::advise_need() {
#if TARGET_OS_MAC
#elif not defined(_MSC_VER)
	posix_fadvise(fileno(file()), 0, buffer_->file_size(), POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED);
#endif
}

size_t decompress(FILE* src, void* dst, size_t dstCapacity, Compressor compressor) {
	switch (compressor) {
	case Compressor::ZLIB:
		return zlib_decompress(src, dst, dstCapacity);
#ifdef WITH_ZSTD
	case Compressor::ZSTD:
		return zstd_decompress(src, dst, dstCapacity);
#endif
	default:
		throw std::runtime_error("Invalid compressor in decompress");
	}
}