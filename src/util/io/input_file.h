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

#ifndef INPUT_FILE_H_
#define INPUT_FILE_H_

#include <vector>
#include <string>
#include <stdexcept>
#include "output_file.h"
#include "deserializer.h"

using std::vector;
using std::string;
using std::runtime_error;

struct InputFile : public Deserializer
{

	enum { BUFFERED = 1 };

	InputFile(const string &file_name, int flags = 0);
	InputFile(OutputFile &tmp_file, int flags = 0);
	void close_and_delete();
	
	string file_name;
	
};

#endif /* BINARY_FILE_H_ */
