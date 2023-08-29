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
#include <string>
#include <stdexcept>
#include <mutex>
#include "stream_entity.h"
#include "exceptions.h"

struct FileSink : public StreamEntity
{
	FileSink(const std::string &file_name, const char *mode = "wb", bool async = false, size_t buffer_size = 0);
#ifndef _MSC_VER
	FileSink(const std::string &file_name, int fd, const char *mode, bool async = false, size_t buffer_size = 0);
#endif
	virtual void close() override;
	virtual void write(const char *ptr, size_t count) override;
	virtual void seek(int64_t p, int origin = SEEK_SET) override;
	virtual void rewind() override;
	virtual int64_t tell() override;
	virtual const std::string& file_name() const override
	{
		return file_name_;
	}
	virtual FILE* file() override
	{
		return f_;
	}
protected:
	FILE *f_;
	const std::string file_name_;
	std::mutex mtx_;
	const bool async_;
	friend struct FileSource;
};