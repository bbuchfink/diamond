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

#ifndef FILE_SINK_H_
#define FILE_SINK_H_

#include <string>
#include <stdexcept>
#include "stream_entity.h"
#include "exceptions.h"

using std::string;

struct FileSink : public StreamEntity
{
	FileSink(const string &file_name, const char *mode = "wb");
#ifndef _MSC_VER
	FileSink(const string &file_name, int fd, const char *mode);
#endif
	virtual void close();
	virtual void write(const char *ptr, size_t count);
	virtual void seek(size_t p);
	virtual void rewind();
	virtual size_t tell();
	virtual const string& file_name() const
	{
		return file_name_;
	}
	virtual FILE* file()
	{
		return f_;
	}
protected:
	FILE *f_;
	const string file_name_;
	friend struct FileSource;
};

#endif