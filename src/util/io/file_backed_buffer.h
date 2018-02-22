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

#ifndef FILE_BACKED_BUFFER_H_
#define FILE_BACKED_BUFFER_H_

#include "serializer.h"
#include "temp_file.h"

struct FileBackedBuffer : public Serializer, public Deserializer
{

	FileBackedBuffer():
		f_(TempFile()),
		Serializer(f_)
	{}

	void rewind()
	{
		Deserializer::f_ = new InputFile(f_, InputFile::BUFFERED);
	}

	~FileBackedBuffer()
	{
		if (Deserializer::f_ == NULL)
			rewind();
		Deserializer::f_->close_and_delete();
		delete Deserializer::f_;
	}

private:

	TempFile f_;

};

#endif
