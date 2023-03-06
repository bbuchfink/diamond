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
#include <utility>
#include "exceptions.h"

struct StreamEntity
{
	StreamEntity(bool seekable = false):
		prev_(NULL),
		seekable_(seekable)
	{}
	StreamEntity(StreamEntity *prev, bool seekable = false):
		prev_(prev),
		seekable_(seekable)
	{}
	virtual void rewind()
	{
		throw UnsupportedOperation();
	}
	virtual void seek(int64_t p, int origin)
	{
		throw UnsupportedOperation();
	}
	virtual int64_t tell()
	{
		throw UnsupportedOperation();
	}
	virtual size_t read(char *ptr, size_t count)
	{
		throw UnsupportedOperation();
	}
	virtual std::pair<const char*, const char*> read()
	{
		throw UnsupportedOperation();
	}
	virtual void close()
	{
		prev_->close();
	}
	virtual const std::string& file_name() const
	{
		return prev_->file_name();
	}
	virtual void write(const char *ptr, size_t count)
	{
		throw UnsupportedOperation();
	}
	virtual std::pair<char*, char*> write_buffer()
	{
		throw UnsupportedOperation();
	}
	virtual void flush(size_t count)
	{
		throw UnsupportedOperation();
	}
	virtual void putback(const char *p, size_t count)
	{
		throw UnsupportedOperation();
	}
	virtual int64_t file_size() {
		if (prev_)
			return prev_->file_size();
		throw UnsupportedOperation();
	}
	virtual FILE* file()
	{
		return prev_->file();
	}
	virtual StreamEntity* root()
	{
		return prev_ ? prev_->root() : this;
	}
	virtual ~StreamEntity()
	{
		delete prev_;
	}
	bool seekable() const {
		return seekable_;
	}
protected:
	StreamEntity *prev_;
	bool seekable_;
};
