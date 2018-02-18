/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef COMPRESSED_STREAM_H_
#define COMPRESSED_STREAM_H_

#include <string>
#include <memory>
#include <zlib.h>
#include "sink.h"
#include "source.h"
#include "../util.h"

using std::string;
using std::auto_ptr;

struct Stream_read_exception : public std::runtime_error
{
	Stream_read_exception(size_t line_count, const char *msg) :
		runtime_error(string("Error reading input stream at line ") + to_string(line_count) + ": " + msg)
	{}
};

struct ZlibSource : public Source
{
	ZlibSource(Source *source);
	virtual size_t read(char *ptr, size_t count);
	virtual void close();
	virtual void rewind();
	virtual const string& file_name() const
	{
		return source_->file_name();
	}
	virtual ~ZlibSource()
	{
		delete source_;
	}
private:
	void init();
	Source *source_;
	z_stream strm;
	static const size_t chunk_size = 1llu << 20;
	auto_ptr<char> in, out;
	size_t read_, total_;
	bool eos_;
};

struct ZlibSink : public Sink
{
	ZlibSink(Sink *sink);
	virtual void close();
	virtual void write(const char *ptr, size_t count);
	virtual FILE* file()
	{
		return sink_->file();
	}
	virtual ~ZlibSink()
	{
		delete sink_;
	}
private:
	Sink *sink_;
	void deflate_loop(const char *ptr, size_t count, int code);
	static const size_t chunk_size = 1llu << 20;
	z_stream strm;
	auto_ptr<char> out;
};

#endif