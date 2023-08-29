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
#include <utility>
#include <memory>
#include "stream_entity.h"

struct OutputStreamBuffer : public StreamEntity
{
	OutputStreamBuffer(StreamEntity* prev);
	virtual std::pair<char*, char*> write_buffer() override;
	virtual void flush(size_t count) override;
	virtual void seek(int64_t p, int origin) override;
	virtual void rewind() override;
	virtual int64_t tell() override;
private:

	static const size_t STDOUT_BUF_SIZE = 4096;
	const size_t buf_size_;
	std::unique_ptr<char[]> buf_;
};