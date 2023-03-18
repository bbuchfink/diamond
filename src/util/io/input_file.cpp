/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include <algorithm>
#include <iostream>
#include "../util/algo/MurmurHash3.h"
#ifndef _MSC_VER
#include <stdio.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include "input_file.h"
#include "file_source.h"
#include "compressed_stream.h"
#include "input_stream_buffer.h"
#include "../../basic/config.h"
#include "temp_file.h"
#ifdef WITH_ZSTD
#include "zstd_stream.h"
#endif

using std::string;

static Compressor detect_compressor(const char* b) {
	if ((b[0] == '\x1F' && b[1] == '\x8B')         // gzip header
		|| (b[0] == '\x78' && (b[1] == '\x01'      // zlib header
			|| b[1] == '\x9C'
			|| b[1] == '\xDA')))
		return Compressor::ZLIB;
	if (b[0] == '\x28' && b[1] == '\xb5' && b[2] == '\x2f' && b[3] == '\xfd')
		return Compressor::ZSTD;
	return Compressor::NONE;
}

static StreamEntity* make_decompressor(const Compressor c, StreamEntity* buffer) {
	switch (c) {
	case Compressor::ZLIB:
		return new ZlibSource(buffer);
	case Compressor::ZSTD:
#ifdef WITH_ZSTD
		return new ZstdSource(buffer);
#else
		throw std::runtime_error("Executable was not compiled with ZStd support.");
#endif
	default:
		throw std::runtime_error("");
	}
}

InputFile::InputFile(const string &file_name, int flags) :
	Deserializer(new InputStreamBuffer(new FileSource(file_name), flags)),
	file_name(file_name),
	unlinked(false),
	temp_file(false)
{
	if (file_name.empty() || file_name == "-")
		return;
#ifndef _MSC_VER
	struct stat buf;
	if (stat(file_name.c_str(), &buf) < 0) {
		perror(0);
		throw std::runtime_error(string("Error calling stat on file ") + file_name);
	}
	if (!S_ISREG(buf.st_mode))
		return;
#endif
	if (flags & NO_AUTODETECT)
		return;
	FileSource *source = dynamic_cast<FileSource*>(buffer_->root());
	char b[4];
	size_t n = source->read(b, 4);
	/*if (n == 2)
		source->putback(b[1]);
	if (n >= 1)
		source->putback(b[0]);*/
	buffer_->putback(b, n);
	if (n < 4)
		return;
	const auto c = detect_compressor(b);
	if(c != Compressor::NONE)
		buffer_ = new InputStreamBuffer(make_decompressor(c, buffer_));
}

InputFile::InputFile(TempFile &tmp_file, int flags) :
	Deserializer(new InputStreamBuffer(new FileSource(tmp_file.file_name(), tmp_file.file()), flags)),
	file_name(tmp_file.file_name()),
	unlinked(tmp_file.unlinked),
	temp_file(true)
{
	tmp_file.rewind();
}

InputFile::InputFile(OutputFile& tmp_file, int flags) :
	Deserializer(new InputStreamBuffer(new FileSource(tmp_file.file_name(), tmp_file.file()), flags)),
	file_name(tmp_file.file_name()),
	unlinked(false),
	temp_file(true)
{
	tmp_file.rewind();
}

void InputFile::close_and_delete()
{
	close();
	if (!unlinked && remove(file_name.c_str()) != 0)
		std::cerr << "Warning: Failed to delete temporary file " << file_name << std::endl;
}

uint64_t InputFile::hash() {
	const size_t SIZE = 4096;
	char buf[SIZE];
	size_t n;
	char h[16];
	std::fill(h, h + 16, '\0');
	while ((n = read_raw(buf, SIZE)) > 0)
		MurmurHash3_x64_128(buf, (int)n, h, h);
	uint64_t r;
	memcpy(&r, h, 8);
	return r;
}
