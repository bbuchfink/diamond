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

#ifndef FILE_BUFFER_H_
#define FILE_BUFFER_H_

#include <string>
#include "../io/temp_file.h"

using std::string;

struct FileBuffer
{

	FileBuffer():
		out_(TempFile()),
		in_(NULL)
	{}

	FileBuffer& operator<<(int x)
	{
		out_.write(&x, 1);
		return *this;
	}

	FileBuffer& operator<<(const string &s)
	{
		out_.write(s.c_str(), s.length() + 1);
		return *this;
	}

	FileBuffer& operator>>(int &x)
	{
		in_->read(&x, 1);
		return *this;
	}

	FileBuffer& operator>>(string &s)
	{
		in_->read_until(s, '\0');
		return *this;
	}

	void rewind()
	{
		in_ = new InputFile(out_, InputFile::BUFFERED);
	}

	~FileBuffer()
	{
		if (!in_)
			rewind();
		in_->close_and_delete();
		delete in_;
	}

private:

	TempFile out_;
	InputFile *in_;

};

#endif