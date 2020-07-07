/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
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
#include <mutex>
#include "temp_file.h"

struct AsyncFile : public TempFile {

	AsyncFile():
		TempFile()
	{}

	template<typename _t>
	void write(const _t *ptr, size_t count)
	{
		std::lock_guard<std::mutex> guard(mtx_);
		write_raw((const char*)ptr, count * sizeof(_t));
	}

private:

	std::mutex mtx_;

};