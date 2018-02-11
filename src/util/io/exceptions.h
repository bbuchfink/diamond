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

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <stdexcept>
#include <string>

using std::string;

struct UnsupportedOperation : public std::runtime_error
{
	UnsupportedOperation() :
		std::runtime_error("Unsupported I/O operation.")
	{ }
};

struct File_open_exception : public std::runtime_error
{
	File_open_exception(const string &file_name) :
		std::runtime_error(string("Error opening file " + file_name))
	{ }
};

struct File_read_exception : public std::runtime_error
{
	File_read_exception(const string &file_name) :
		runtime_error(string("Error reading file ") + file_name)
	{ }
};

struct File_write_exception : public std::runtime_error
{
	File_write_exception(const string &file_name) :
		runtime_error(string("Error writing file ") + file_name)
	{ }
};

#endif