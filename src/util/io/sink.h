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

#ifndef SINK_H_
#define SINK_H_

#include "exceptions.h"

struct Sink
{
	virtual void close() = 0;
	virtual void write(const char *ptr, size_t count) = 0;
	virtual void seekp(size_t p)
	{
		throw UnsupportedOperation();
	}
	virtual size_t tell()
	{
		throw UnsupportedOperation();
	}
	virtual FILE* file() = 0;
	virtual ~Sink()
	{}
};

#endif