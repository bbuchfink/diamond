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

#include "text_input_file.h"

using std::string;

TextInputFile::TextInputFile(const string &file_name, char line_separator) :
	InputFile(file_name),
	line_count(0),
	putback_line_(false),
	eof_(false),
	line_separator(line_separator)
{
}

TextInputFile::TextInputFile(TempFile &tmp_file, char line_separator) :
	InputFile(tmp_file),
	line_count(0),
	putback_line_(false),
	eof_(false),
	line_separator(line_separator)
{}

TextInputFile::TextInputFile(OutputFile& out_file, char line_separator) :
	InputFile(out_file),
	line_count(0),
	putback_line_(false),
	eof_(false),
	line_separator(line_separator)
{}

void TextInputFile::rewind()
{
	InputFile::rewind();
	line_count = 0;
	putback_line_ = false;
	eof_ = false;
	line.clear();
}

bool TextInputFile::eof() const
{
	return eof_;
}

void TextInputFile::getline()
{
	//assert(strlen(line_separator) <= 2);
	if (putback_line_) {
		putback_line_ = false;
		++line_count;
		return;
	}
	line.clear();
	//if (strlen(line_separator) == 1) {
	eof_ = !read_to(std::back_inserter(line), line_separator);
	++line_count;
	/* }
	else {
		auto r = read_to(std::back_inserter(line), line_separator[0], line_separator[1]);
		eof_ = !r.first;
		line_count += r.second;
	}*/
	if (!line.empty() && line.back() == '\r')
		line.resize(line.length() - 1);
}

void TextInputFile::putback_line()
{
	putback_line_ = true;
	--line_count;
}

/*void TextInputFile::set_separator(char separator) {
	this->line_separator = separator;
}*/