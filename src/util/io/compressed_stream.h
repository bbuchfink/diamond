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

#pragma once
#include <string>
#include <zlib.h>
#include "stream_entity.h"

struct ZlibSource : public StreamEntity
{
	ZlibSource(StreamEntity *prev);
	virtual size_t read(char *ptr, size_t count);
	virtual void close();
	virtual void rewind();
private:
	void init();
	z_stream strm;
	static const size_t chunk_size = 1llu << 20;
	bool eos_;
};

struct ZlibSink : public StreamEntity
{
	ZlibSink(StreamEntity *prev);
	virtual void close();
	virtual void write(const char *ptr, size_t count);
private:
	void deflate_loop(const char *ptr, size_t count, int code);
	static const size_t chunk_size = 1llu << 20;
	z_stream strm;
};