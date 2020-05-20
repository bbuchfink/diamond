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

#ifndef FILE_SOURCE_H_
#define FILE_SOURCE_H_

#include "stream_entity.h"

struct FileSource : public StreamEntity
{
	FileSource(const string &file_name);
	FileSource(const string &file_name, FILE *file);
	virtual void rewind();
	virtual void seek(size_t pos);
	virtual void seek_forward(size_t n);
	virtual size_t read(char *ptr, size_t count);
	virtual void close();
	virtual const string& file_name() const
	{
		return file_name_;
	}
	virtual FILE* file()
	{
		return f_;
	}
	//void putback(char c);
	~FileSource()
	{}
protected:
	FILE *f_;
	const string file_name_;
};

#endif