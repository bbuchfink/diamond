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

#pragma once
#include <utility>
#include <memory>
#include <thread>
#include "stream_entity.h"

struct InputStreamBuffer : public StreamEntity
{
	enum { ASYNC = 4 };
	InputStreamBuffer(StreamEntity* prev, int flags = 0);
	virtual void rewind() override;
	virtual void seek(int64_t p, int origin) override;
	virtual std::pair<const char*, const char*> read() override;
	virtual void putback(const char* p, size_t n) override;
	virtual void close() override;
	virtual int64_t tell() override;
private:

	static void load_worker(InputStreamBuffer *buf);
	
	const size_t buf_size_;
	std::unique_ptr<char[]> buf_, load_buf_;
	size_t putback_count_, load_count_, file_offset_;
	bool async_;
	std::thread* load_worker_;
};