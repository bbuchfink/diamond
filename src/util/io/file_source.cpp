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

#include <iostream>
#include <stdio.h>
#ifdef _MSC_VER
#define NOMINMAX
#include <Windows.h>
#else
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "file_source.h"
#include "../system.h"

using std::endl;
using std::string;
using std::runtime_error;

FileSource::FileSource(const string &file_name) :
	StreamEntity(true),
	file_name_(file_name)
{
	static const char* msg = "\nError opening file ";
	const bool is_stdin = file_name.empty() || file_name == "-";
#ifdef _MSC_VER
	f_ = is_stdin ? stdin : fopen(file_name.c_str(), "rb");
#else

	struct stat buf;
	if (!is_stdin && stat(file_name.c_str(), &buf) < 0) {
		perror((msg + file_name).c_str());
		throw FileOpenException(file_name);
	}
	if (is_stdin || !S_ISREG(buf.st_mode))
		seekable_ = false;

	int fd_ = is_stdin ? 0 : POSIX_OPEN2(file_name.c_str(), O_RDONLY);
	if (fd_ < 0) {
		perror((msg + file_name).c_str());
		throw FileOpenException(file_name);
	}
	f_ = fdopen(fd_, "rb");
#endif
	if (f_ == 0) {
		perror((msg + file_name).c_str());
		throw FileOpenException(file_name);
	}
}

FileSource::FileSource(const string &file_name, FILE *file):
	StreamEntity(false),
	f_(file),
	file_name_(file_name)
{
}

void FileSource::rewind()
{
	::rewind(f_);
}

void FileSource::seek(int64_t pos, int origin)
{
#ifdef _MSC_VER
	if (_fseeki64(f_, pos, SEEK_SET) != 0) {
		perror(0);
		throw std::runtime_error("Error executing seek on file " + file_name_);
	}
#else
	if (fseek(f_, pos, SEEK_SET) < 0) {
		perror(0);
		throw std::runtime_error("Error calling fseek.");
	}
#endif
}

void FileSource::seek_forward(size_t n)
{
#ifdef _MSC_VER
	if (_fseeki64(f_, (int64_t)n, SEEK_CUR) != 0) {
		perror(0);
		throw std::runtime_error("Error executing seek on file " + file_name_);
	}
#else
	if (fseek(f_, n, SEEK_CUR) < 0) {
		perror(0);
		throw std::runtime_error("Error calling fseek.");
	}
#endif
}

/*void FileSource::putback(char c)
{
	if ((char)ungetc((int)c, f_) != c) {
		perror(0);
		throw std::runtime_error("Error calling ungetc.");
	}
}*/

size_t FileSource::read(char *ptr, size_t count)
{
	size_t n;
	if ((n = fread(ptr, 1, count, f_)) != count) {
		if (feof(f_) != 0)
			return n;
		else {
			perror(0);
			throw File_read_exception(file_name_);
		}
	}
	return n;
}

void FileSource::close()
{
	if (f_) {
		if (fclose(f_) != 0) {
			perror(0);
			throw std::runtime_error(string("Error closing file ") + file_name_);
		}
		f_ = 0;
	}
}

int64_t FileSource::tell()
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