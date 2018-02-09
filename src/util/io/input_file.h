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

#ifndef INPUT_FILE_H_
#define INPUT_FILE_H_

#include <vector>
#include <string>
#include <stdexcept>
#include "output_file.h"
#include "source.h"

using std::vector;
using std::string;
using std::runtime_error;

struct InputFile
{
	InputFile(const string &file_name);
	void rewind();
	InputFile(const OutputFile &tmp_file);
	void seek(size_t pos);
	void seek_forward(size_t n);
	template<class _t>
	size_t read(_t *ptr, size_t count)
	{
		return source_->read((char*)ptr, count*sizeof(_t)) / sizeof(_t);
	}
	void read_c_str(string &s);
	void close_and_delete();
	void close();
	~InputFile();
	
	string file_name;
	
private:
	
	Source *source_;
	
};

#endif /* BINARY_FILE_H_ */
