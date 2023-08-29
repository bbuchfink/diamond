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

#include <stdexcept>
#include "compressed_stream.h"

using std::pair;

void ZlibSource::init()
{
	eos_ = false;
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	int ret = inflateInit2(&strm, 15 + 32);
	if (ret != Z_OK)
		throw std::runtime_error("Error opening compressed file (inflateInit): " + file_name());
}

ZlibSource::ZlibSource(StreamEntity *prev):
	StreamEntity(prev)
{
	init();
}

size_t ZlibSource::read(char *ptr, size_t count)
{
	strm.avail_out = (uInt)count;
	strm.next_out = (Bytef*)ptr;
	while (strm.avail_out > 0 && !eos_) {
		if (strm.avail_in == 0) {
			pair<const char*, const char*> in = prev_->read();

			strm.avail_in = (uInt)(in.second - in.first);
			if (strm.avail_in == 0) {
				eos_ = true;
				break;
			}
			strm.next_in = (Bytef*)in.first;
		}

		int ret = inflate(&strm, Z_NO_FLUSH);
		if (ret == Z_STREAM_END) {
			int ret = inflateInit2(&strm, 15 + 32);
			if (ret != Z_OK)
				throw std::runtime_error("Error initializing compressed stream (inflateInit): " + file_name());
		}
		else if (ret != Z_OK)
			throw std::runtime_error("Inflate error.");
	}
	return count - strm.avail_out;
}

void ZlibSource::close()
{
	inflateEnd(&strm);
	prev_->close();
}

void ZlibSource::rewind()
{
	prev_->rewind();
	inflateEnd(&strm);
	init();
}

ZlibSink::ZlibSink(StreamEntity *prev):
	StreamEntity(prev)
{
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK)
		throw std::runtime_error("deflateInit error");
}

void ZlibSink::deflate_loop(const char *ptr, size_t count, int flush)
{
	strm.avail_in = (uInt)count;
	strm.next_in = (Bytef*)ptr;
	do {
		pair<char*, char*> out = prev_->write_buffer();
		const size_t chunk_size = out.second - out.first;
		strm.avail_out = (uInt)chunk_size;
		strm.next_out = (Bytef*)out.first;
		if (deflate(&strm, flush) == Z_STREAM_ERROR)
			throw std::runtime_error("deflate error");
		prev_->flush(chunk_size - strm.avail_out);
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
	prev_->close();
}