/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include <stdexcept>
#include "compressed_stream.h"

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

Compressed_istream::Compressed_istream(const string & file_name) :
	file_name_(file_name),
	line_count(0),
	s_(file_name, std::ios_base::in | std::ios_base::binary),	
	putback_line_(false)
{ }

size_t Compressed_istream::read(char * ptr, size_t count)
{
	s_.read(ptr, count);
	const size_t n = s_.gcount();
	if (n != count) {
		if (s_.eof())
			return n;
		else
			throw std::runtime_error("Error reading file " + file_name_);
	}
	return n;
}

void Compressed_istream::putback(char c)
{
	s_.putback(c);
	if (!s_.good())
		throw std::runtime_error("Error reading file " + file_name_);
}

void Compressed_istream::getline()
{
	if (!putback_line_) {
		std::getline(s_, line);
		if (!s_.good() && !s_.eof())
			throw Stream_read_exception(line_count, "I/O error");
		const size_t s = line.length() - 1;
		if (!line.empty() && line[s] == '\r')
			line.resize(s);
	}
	else
		putback_line_ = false;
	++line_count;
}

void Compressed_istream::putback_line()
{
	putback_line_ = true;
	--line_count;
}

Compressed_ostream::Compressed_ostream(const string &file_name):
	Output_stream(file_name),
	out(new char[chunk_size])
{
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK)
		throw std::runtime_error("deflateInit error");
}

void Compressed_ostream::deflate_loop(const char * ptr, size_t count, int flush)
{
	strm.avail_in = (uInt)count;
	strm.next_in = (Bytef*)ptr;
	do {
		strm.avail_out = chunk_size;
		strm.next_out = (Bytef*)out.get();
		if (deflate(&strm, flush) == Z_STREAM_ERROR)
			throw std::runtime_error("deflate error");
		const size_t have = chunk_size - strm.avail_out;
		Output_stream::write(out.get(), have);
	} while (strm.avail_out == 0);
}

void Compressed_ostream::write(const char * ptr, size_t count)
{
	deflate_loop(ptr, count, Z_NO_FLUSH);
}

void Compressed_ostream::close()
{
	deflate_loop(0, 0, Z_FINISH);
	deflateEnd(&strm);
}