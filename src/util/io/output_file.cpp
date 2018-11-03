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

#include <iostream>
#include <stdio.h>
#include "output_file.h"
#include "file_sink.h"
#include "output_stream_buffer.h"
#include "compressed_stream.h"

OutputFile::OutputFile(const string &file_name, bool compressed, const char *mode) :
	Serializer(new OutputStreamBuffer(new FileSink(file_name, mode))),
	file_name_(file_name)
{
	if (compressed) {
		buffer_ = new OutputStreamBuffer(new ZlibSink(buffer_));
		reset_buffer();
	}
}

#ifndef _MSC_VER
OutputFile::OutputFile(pair<string, int> fd, const char *mode):
	Serializer(new OutputStreamBuffer(new FileSink(fd.first, fd.second, mode))),
	file_name_(fd.first)
{
}
#endif

void OutputFile::remove()
{
	if (::remove(file_name_.c_str()) != 0)
		std::cerr << "Warning: Failed to delete file " << file_name_ << std::endl;
}