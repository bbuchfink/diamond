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

#include <iostream>
#include "../lib/murmurhash/MurmurHash3.h"
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
#include "temp_file.h"
#ifdef WITH_ZSTD
#include "zstd_stream.h"
#endif

using std::string;
using std::runtime_error;

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

static StreamEntity* make_decompressor(const Compressor c, InputStreamBuffer* buffer) {
	switch (c) {
	case Compressor::ZLIB:
		return new ZlibSource(buffer);
	case Compressor::ZSTD:
#ifdef WITH_ZSTD
		return new ZstdSource(buffer);
#else
		throw runtime_error("Executable was not compiled with ZStd support.");
#endif
	default:
		throw runtime_error("");
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
	const string begin = peek(4);
	if (begin.length() < 4)
		return;
	const Compressor c = detect_compressor(begin.c_str());
	if (c != Compressor::NONE)
		buffer_ = new InputStreamBuffer(make_decompressor(c, buffer_));
}

InputFile::InputFile(TempFile &tmp_file, int flags, Compressor c) :
	Deserializer(new InputStreamBuffer(new FileSource(tmp_file.file_name(), tmp_file.file()), flags)),
	file_name(tmp_file.file_name()),
	unlinked(tmp_file.unlinked),
	temp_file(true)
{
	tmp_file.rewind();
	if (c != Compressor::NONE)
		buffer_ = new InputStreamBuffer(make_decompressor(c, buffer_));
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