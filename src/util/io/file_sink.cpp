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

#include <iostream>
#include <stdio.h>
#ifdef _MSC_VER
#define NOMINMAX
#include <Windows.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#endif

#include "file_sink.h"
#include "../system.h"

using std::endl;
using std::string;
using std::runtime_error;

FileSink::FileSink(const string &file_name, const char *mode):
	file_name_(file_name)
{
#ifdef _MSC_VER
	f_ = file_name.length() == 0 ? stdout : fopen(file_name.c_str(), mode);
#else
	int fd_ = file_name.length() == 0 ? 1 : POSIX_OPEN(file_name.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
	if (fd_ < 0) {
		perror(0);
		throw File_open_exception(file_name_);
	}
	f_ = fdopen(fd_, mode);
#endif	
	if (f_ == 0) {
		perror(0);
		throw File_open_exception(file_name);
	}
}

#ifndef _MSC_VER
FileSink::FileSink(const string &file_name, int fd, const char *mode):
	f_(fdopen(fd, mode)),
	file_name_(file_name)	
{
	if (f_ == 0) {
		perror(0);
		throw File_open_exception(file_name);
	}
}
#endif

void FileSink::close()
{
	if (f_ && f_ != stdout) {
		if (fclose(f_) != 0) {
			perror(0);
			throw std::runtime_error(string("Error closing file ") + file_name_);
		}
		f_ = 0;
	}
}

void FileSink::write(const char *ptr, size_t count)
{
	size_t n;
	if ((n = fwrite((const void*)ptr, 1, count, f_)) != count) {
		perror(0);
		throw File_write_exception(file_name_);
	}
}

void FileSink::seekp(size_t p)
{
#ifdef _MSC_VER
	if (_fseeki64(f_, (int64_t)p, SEEK_SET) != 0) {
		perror(0);
		throw std::runtime_error("Error calling fseek.");
	}
#else
	if (fseek(f_, p, SEEK_SET) != 0) {
		perror(0);
		throw std::runtime_error("Error calling fseek.");
	}
#endif
}

size_t FileSink::tell()
{
#ifdef _MSC_VER
	int64_t x;
	if ((x = _ftelli64(f_)) == (int64_t)-1)
		throw std::runtime_error("Error executing ftell on stream " + file_name_);
	return (size_t)x;
#else
	const long n = ftell(f_);
	if (n < 0) {
		perror(0);
		throw std::runtime_error("Error calling ftell.");
	}
	return n;
#endif
}