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

#include <stdexcept>
#include <stdio.h>
#ifndef _MSC_VER
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include "compressed_stream.h"

ZlibSource::ZlibSource(Source *source):
	source_(source),
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
		throw std::runtime_error("Error opening compressed file (inflateInit): " + file_name());
}

size_t ZlibSource::read(char *ptr, size_t count)
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
			strm.avail_in = (uInt)source_->read(in.get(), chunk_size);
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
			int ret = inflateInit2(&strm, 15 + 32);
			if (ret != Z_OK)
				throw std::runtime_error("Error initializing compressed stream (inflateInit): " + file_name());
		}
		else if (ret != Z_OK)
			throw std::runtime_error("Inflate error.");

		read_ = 0;
		total_ = chunk_size - strm.avail_out;
	} while (true);
}

void ZlibSource::close()
{
	inflateEnd(&strm);
	source_->close();
}

ZlibSink::ZlibSink(Sink* sink):
	sink_(sink),
	out(new char[chunk_size])
{
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK)
		throw std::runtime_error("deflateInit error");
}

void ZlibSink::deflate_loop(const char * ptr, size_t count, int flush)
{
	strm.avail_in = (uInt)count;
	strm.next_in = (Bytef*)ptr;
	do {
		strm.avail_out = chunk_size;
		strm.next_out = (Bytef*)out.get();
		if (deflate(&strm, flush) == Z_STREAM_ERROR)
			throw std::runtime_error("deflate error");
		const size_t have = chunk_size - strm.avail_out;
		sink_->write(out.get(), have);
	} while (strm.avail_out == 0);
}

void ZlibSink::write(const char * ptr, size_t count)
{
	deflate_loop(ptr, count, Z_NO_FLUSH);
}

void ZlibSink::close()
{
	deflate_loop(0, 0, Z_FINISH);
	deflateEnd(&strm);
	sink_->close();
}