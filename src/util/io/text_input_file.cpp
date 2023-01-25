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

TextInputFile::TextInputFile(const string &file_name) :
	InputFile(file_name),
	line_count(0),
	line_buf_used_(0),
	line_buf_end_(0),
	putback_line_(false),
	eof_(false)
{
}

TextInputFile::TextInputFile(TempFile &tmp_file) :
	InputFile(tmp_file),
	line_count(0),
	line_buf_used_(0),
	line_buf_end_(0),
	putback_line_(false),
	eof_(false)
{}

TextInputFile::TextInputFile(OutputFile& out_file) :
	InputFile(out_file),
	line_count(0),
	line_buf_used_(0),
	line_buf_end_(0),
	putback_line_(false),
	eof_(false)
{}

void TextInputFile::rewind()
{
	InputFile::rewind();
	line_count = 0;
	line_buf_used_ = 0;
	line_buf_end_ = 0;
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
	if (putback_line_) {
		putback_line_ = false;
		++line_count;
		return;
	}
	line.clear();
	while (true) {
		const char *p = (const char*)memchr(&line_buf_[line_buf_used_], '\n', line_buf_end_ - line_buf_used_);
		if (p == 0) {
			line.append(&line_buf_[line_buf_used_], line_buf_end_ - line_buf_used_);
			line_buf_end_ = read(line_buf_, line_buf_size);
			line_buf_used_ = 0;
			if (line_buf_end_ == 0) {
				eof_ = true;
				++line_count;
				return;
			}
		}
		else {
			const size_t n = (p - line_buf_) - line_buf_used_;
			line.append(&line_buf_[line_buf_used_], n);
			line_buf_used_ += n + 1;
			const size_t s = line.length() - 1;
			if (!line.empty() && line[s] == '\r')
				line.resize(s);
			++line_count;
			return;
		}
	}
}

void TextInputFile::putback_line()
{
	putback_line_ = true;
	--line_count;
}
