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
#include <stdexcept>
#include <string>

struct UnsupportedOperation : public std::runtime_error
{
	UnsupportedOperation() :
		std::runtime_error("Unsupported I/O operation.")
	{ }
};

struct FileOpenException : public std::runtime_error
{
	FileOpenException(const std::string &file_name) :
		std::runtime_error(std::string("Error opening file " + file_name))
	{ }
};

struct File_read_exception : public std::runtime_error
{
	File_read_exception(const std::string &file_name) :
		runtime_error(std::string("Error reading file ") + file_name)
	{ }
};

struct File_write_exception : public std::runtime_error
{
	File_write_exception(const std::string &file_name) :
		runtime_error(std::string("Error writing file ") + file_name)
	{ }
};

struct EndOfStream : public std::runtime_error
{
	EndOfStream() :
		runtime_error("Unexpected end of input.")
	{ }
};

struct StreamReadException : public std::runtime_error
{
	StreamReadException(size_t line_count, const char *msg) :
		runtime_error(std::string("Error reading input stream at line ") + std::to_string(line_count) + ": " + msg)
	{}
};