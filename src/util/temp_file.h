/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef TEMP_FILE_H_
#define TEMP_FILE_H_

#include "binary_file.h"

struct Temp_file : public Output_stream
{

	Temp_file();
#ifndef _MSC_VER
	virtual ~Temp_file()
	{}
#endif
	static string get_temp_dir();
	static unsigned n;
	static uint64_t hash_key;

};

#endif /* TEMP_FILE_H_ */
