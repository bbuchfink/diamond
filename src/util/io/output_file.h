/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <string>
#include "serializer.h"

enum class Compressor { NONE, ZLIB, ZSTD };

struct TempFileData {
	std::string name;
	int fd;
	bool unlinked;
};

struct OutputFile : public Serializer
{
	OutputFile(const std::string &file_name, Compressor compressor = Compressor::NONE, const char *mode = "wb");
	OutputFile(const TempFileData& d, Compressor compressor, const char *mode);

	void remove();

	std::string file_name() const
	{
		return file_name_;
	}

	void advise_need();

protected:

	std::string file_name_;

};

size_t decompress(FILE* src, void* dst, size_t dstCapacity, Compressor compressor);