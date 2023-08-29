/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include <vector>
#include <string>
#include <stdexcept>
#include <stdint.h>
#include "output_file.h"
#include "deserializer.h"
#include "temp_file.h"

const size_t MEGABYTES = 1 << 20;
const size_t GIGABYTES = 1 << 30;
const size_t KILOBYTES = 1 << 10;

struct InputFile : public Deserializer
{

	enum { BUFFERED = 1, NO_AUTODETECT = 2 };

	InputFile(const std::string &file_name, int flags = 0);
	InputFile(TempFile &tmp_file, int flags = 0);
	InputFile(OutputFile& tmp_file, int flags = 0);
	void close_and_delete();
	uint64_t hash();
	
	std::string file_name;
	bool unlinked, temp_file;
	
};
