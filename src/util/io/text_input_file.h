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

#pragma once
#include <string>
#include "input_file.h"

struct TextInputFile : public InputFile
{
	/*TextInputFile(const std::string& file_name, const char* line_separator = "\n");
	TextInputFile(TempFile &tmp_file, const char* line_separator = "\n");
	TextInputFile(OutputFile& out_file, const char* line_separator = "\n");*/
	TextInputFile(const std::string& file_name, char line_separator = '\n');
	TextInputFile(TempFile& tmp_file, char line_separator = '\n');
	TextInputFile(OutputFile& out_file, char line_separator = '\n');
	void rewind();
	bool eof() const;
	void putback(char c);
	void getline();
	void putback_line();
	//void set_separator(const char* separator);
	operator bool() const {
		return !eof();
	}

	std::string line;
	size_t line_count;

protected:

	bool putback_line_, eof_;
	const char line_separator;

};
