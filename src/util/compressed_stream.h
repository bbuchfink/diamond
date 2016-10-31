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

#ifndef COMPRESSED_STREAM_H_
#define COMPRESSED_STREAM_H_

#include <string>
#include <memory>
#include <zlib.h>
#include "binary_file.h"
#include "util.h"

using std::string;
using std::auto_ptr;

struct Stream_read_exception : public std::runtime_error
{
	Stream_read_exception(size_t line_count, const char *msg) :
		runtime_error(string("Error reading input stream at line ") + to_string(line_count) + ": " + msg)
	{}
};

struct Compressed_istream : public Input_stream
{
	Compressed_istream(const string &file_name);
	virtual size_t read_bytes(char *ptr, size_t count);
	virtual void close();
	static Input_stream* auto_detect(const string &file_name);
private:
	z_stream strm;
	static const size_t chunk_size = 1llu << 20;
	auto_ptr<char> in, out;
	size_t read_, total_;
	bool eos_;
};

struct Compressed_ostream : public Output_stream
{
	Compressed_ostream(const string &file_name);
#ifndef _MSC_VER
	virtual ~Compressed_ostream()
	{}
#endif
	virtual void write(const char *ptr, size_t count);
	virtual void close();
private:
	void deflate_loop(const char *ptr, size_t count, int code);
	static const size_t chunk_size = 1llu << 20;
	z_stream strm;
	auto_ptr<char> out;
};

#endif