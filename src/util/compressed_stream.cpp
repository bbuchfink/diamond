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

Compressed_istream::Compressed_istream(const string &file_name):
	Input_stream(file_name),
	in(new char[chunk_size]),
	out(new char[chunk_size]),
	read_(0),
	total_(0),
	eos_(false)
{
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	strm.avail_out = chunk_size;
	strm.next_out = (Bytef*)out.get();
	int ret = inflateInit2(&strm, 15 + 32);
	if (ret != Z_OK)
		throw std::runtime_error("Error opening compressed file (inflateInit): " + file_name);
}

size_t Compressed_istream::read_bytes(char *ptr, size_t count)
{
	size_t n = 0;
	do {
		size_t m = std::min(count - n, total_ - read_);
		memcpy(ptr, &out.get()[read_], m);
		read_ += m;
		ptr += m;
		n += m;
		if (count == n || eos_)
			return n;

		if (strm.avail_out > 0 && strm.avail_in == 0) {
			strm.avail_in = (uInt)Input_stream::read_bytes(in.get(), chunk_size);
			if (strm.avail_in == 0) {
				eos_ = true;
				return n;
			}
			strm.next_in = (Bytef*)in.get();
		}

		strm.avail_out = chunk_size;
		strm.next_out = (Bytef*)out.get();

		int ret = inflate(&strm, Z_NO_FLUSH);
		if (ret == Z_STREAM_END) {
			if(strm.avail_in == 0 && feof(this->f_))
				eos_ = true;
			else {
				int ret = inflateInit2(&strm, 15 + 32);
				if (ret != Z_OK)
					throw std::runtime_error("Error initializing compressed stream (inflateInit): " + file_name);
			}
		}
		else if (ret != Z_OK)
			throw std::runtime_error("Inflate error.");

		read_ = 0;
		total_ = chunk_size - strm.avail_out;
	} while (true);
}

void Compressed_istream::close()
{
	inflateEnd(&strm);
	Input_stream::close();
}

bool is_gzip_stream(const unsigned char *b)
{
	return (b[0] == 0x1F && b[1] == 0x8B)         // gzip header
		|| (b[0] == 0x78 && (b[1] == 0x01      // zlib header
			|| b[1] == 0x9C
			|| b[1] == 0xDA));
}

Input_stream *Compressed_istream::auto_detect(const string &file_name)
{
	unsigned char b[2];
	Input_stream f(file_name);
	size_t n = f.read(b, 2);
	f.close();
	if (n == 2 && is_gzip_stream(b))
		return new Compressed_istream(file_name);
	else
		return new Input_stream(file_name);
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