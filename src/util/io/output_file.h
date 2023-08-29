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
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include "serializer.h"
#include "../text_buffer.h"

enum class Compressor { NONE, ZLIB, ZSTD };

struct OutputFile : public Serializer
{
	OutputFile(const std::string &file_name, Compressor compressor = Compressor::NONE, const char *mode = "wb");
#ifndef _MSC_VER
	OutputFile(std::pair<std::string, int> fd, const char *mode);
#endif

	void remove();

	template<typename _k, typename _v>
	void write_map_csv(typename std::map<_k, _v>::const_iterator begin, typename std::map<_k, _v>::const_iterator end)
	{
		TextBuffer buf;
		for (typename std::map<_k, _v>::const_iterator i = begin; i != end; ++i) {
			buf << i->first << '\t' << i->second << '\n';
			write(buf.data(), buf.size());
			buf.clear();
		}
	}

	std::string file_name() const
	{
		return file_name_;
	}

protected:

	std::string file_name_;

};
