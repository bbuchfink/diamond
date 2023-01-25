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
#include <string.h>
#endif

#include "file_sink.h"
#include "../system.h"

using std::endl;
using std::string;
using std::runtime_error;

static const char* msg = "\nError opening file ";

#ifndef _MSC_VER
static int posix_flags(const char* mode) {
	if (strcmp(mode, "wb") == 0)
		return O_WRONLY | O_CREAT | O_TRUNC;
	else if (strcmp(mode, "r+b") == 0)
		return O_RDWR;
	else if (strcmp(mode, "w+b") == 0)
		return O_RDWR | O_CREAT | O_TRUNC;
	throw std::runtime_error("Invalid fopen mode.");
}
#endif

FileSink::FileSink(const string &file_name, const char *mode, bool async, size_t buffer_size):
	file_name_(file_name),
	async_(async)
{
#ifdef _MSC_VER
	f_ = file_name.length() == 0 ? stdout : fopen(file_name.c_str(), mode);
#else
	int fd_ = file_name.length() == 0 ? 1 : POSIX_OPEN(file_name.c_str(), posix_flags(mode), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
	if (fd_ < 0) {
		perror((msg + file_name).c_str());
		throw FileOpenException(file_name_);
	}
	f_ = fdopen(fd_, mode);
#endif	
	if (f_ == 0) {
		perror((msg + file_name).c_str());
		throw FileOpenException(file_name);
	}
	if (buffer_size != 0)
		if (setvbuf(f_, nullptr, _IOFBF, buffer_size) != 0)
			throw std::runtime_error("Error calling setvbuf.");
}

#ifndef _MSC_VER
FileSink::FileSink(const string &file_name, int fd, const char *mode, bool async, size_t buffer_size):
	f_(fdopen(fd, mode)),
	file_name_(file_name),
	async_(async)
{
	if (f_ == 0) {
		perror((msg + file_name).c_str());
		throw FileOpenException(file_name);
	}
	if (buffer_size != 0)
		if (setvbuf(f_, nullptr, _IOFBF, buffer_size) != 0)
			throw std::runtime_error("Error calling setvbuf.");
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
	if (async_) mtx_.lock();
	size_t n;
	if ((n = fwrite((const void*)ptr, 1, count, f_)) != count) {
		if (async_) mtx_.unlock();
		perror(0);
		throw File_write_exception(file_name_);
	}
	if (async_) mtx_.unlock();
}

void FileSink::seek(int64_t p, int origin)
{
#ifdef _MSC_VER
	if (_fseeki64(f_, p, origin) != 0) {
		perror(0);
		throw std::runtime_error("Error calling fseek.");
	}
#else
	if (fseek(f_, p, origin) != 0) {
		perror(0);
		throw std::runtime_error("Error calling fseek.");
	}
#endif
}

int64_t FileSink::tell()
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

void FileSink::rewind()
{
	::rewind(f_);
}