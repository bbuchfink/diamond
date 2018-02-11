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

#include <iostream>
#include <stdio.h>
#include "output_file.h"

void OutputFile::remove()
{
	if (::remove(file_name_.c_str()) != 0)
		std::cerr << "Warning: Failed to delete file " << file_name_ << std::endl;
}

void OutputFile::write_c_str(const string &s)
{
	write(s.c_str(), s.length() + 1);
}

void OutputFile::seekp(size_t p)
{
	sink_->seekp(p);
}

size_t OutputFile::tell()
{
	return sink_->tell();
}