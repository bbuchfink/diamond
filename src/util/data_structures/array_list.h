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

#ifndef ARRAY_LIST_H_
#define ARRAY_LIST_H_

#include <vector>
#include "file_buffer.h"

using std::vector;

struct ArrayList : private FileBuffer
{

	ArrayList():
		FileBuffer(),
		entries_(0)
	{}

	void push_back(const vector<string> &v)
	{
		(*this) << (int)v.size();
		for (vector<string>::const_iterator i = v.begin(); i < v.end(); ++i)
			(*this) << *i;
		entries_ += v.size();
	}

	void rewind()
	{
		FileBuffer::rewind();
	}

	size_t entries() const
	{
		return entries_;
	}

private:

	size_t entries_;

};

#endif