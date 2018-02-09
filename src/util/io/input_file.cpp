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
#include "input_file.h"
#include "file_source.h"
#include "compressed_stream.h"

bool is_gzip_stream(const unsigned char *b)
{
	return (b[0] == 0x1F && b[1] == 0x8B)         // gzip header
		|| (b[0] == 0x78 && (b[1] == 0x01      // zlib header
			|| b[1] == 0x9C
			|| b[1] == 0xDA));
}

InputFile::InputFile(const string &file_name) :
	file_name(file_name),
	source_ (new FileSource(file_name))
{
	if (file_name.empty())
		return;
#ifndef _MSC_VER
	struct stat buf;
	if (stat(file_name.c_str(), &buf) < 0) {
		perror(0);
		throw std::runtime_error(string("Error calling stat on file ") + file_name);
	}
	if (!S_ISREG(buf.st_mode))
		return;
#endif	
	char b[2];
	size_t n = source_->read(b, 2);
	if (n == 2)
		source_->putback(b[1]);
	if (n >= 1)
		source_->putback(b[0]);
	if (n == 2 && is_gzip_stream((const unsigned char*)b))
		source_ = new ZlibSource(source_);
}

InputFile::~InputFile()
{
	delete source_;
}

void InputFile::rewind()
{
	source_->rewind();
}

InputFile::InputFile(const OutputFile &tmp_file) :
	file_name(tmp_file.file_name()),
	source_(new FileSource(file_name, tmp_file.file()))
{
	rewind();
}

void InputFile::seek(size_t pos)
{
	source_->seek(pos);
}

void InputFile::seek_forward(size_t n)
{
	source_->seek_forward(n);
}

void InputFile::read_c_str(string &s)
{
	char c;
	s.clear();
	while (true) {
		if (read(&c, 1) != 1)
			throw std::runtime_error("Unexpected end of file.");
		if (c == 0)
			break;
		s += (char)c;
	}
}

void InputFile::close_and_delete()
{
	source_->close();
#ifdef _MSC_VER
	if (remove(file_name.c_str()) != 0)
		std::cerr << "Warning: Failed to delete temporary file " << file_name << std::endl;
#endif
}

void InputFile::close()
{
	source_->close();
}