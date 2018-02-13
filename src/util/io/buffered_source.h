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

#ifndef BUFFERED_SOURCE_H_
#define BUFFERED_SOURCE_H_

#include <assert.h>
#include <vector>
#include <string>
#include "source.h"

using std::vector;
using std::string;

struct BufferedSource : public Source
{
	BufferedSource(Source* source);
	virtual void rewind();
	virtual void seek(size_t pos);
	virtual void seek_forward(size_t n);
	virtual size_t read(char *ptr, size_t count);
	virtual void close();
	virtual const string& file_name() const;
	virtual bool read_until(string &dst, char delimiter);
	virtual bool read_until(vector<char> &dst, char delimiter);
	virtual ~BufferedSource()
	{
		delete source_;
	}
private:
	char* next()
	{
		return &buf_[start_];
	}
	void pop(char *dst, size_t n);
	size_t fetch()
	{
		start_ = 0;
		return avail_ = source_->read(buf_, BUF_SIZE);
	}

	enum { BUF_SIZE = 4096 };

	Source *source_;
	char buf_[BUF_SIZE];
	size_t start_, avail_;
	bool eof_;
};

#endif