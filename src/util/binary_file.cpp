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

#include <sstream>
#ifdef _MSC_VER
#else
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#endif

#include "../basic/config.h"
#include "binary_file.h"
#include "temp_file.h"
#include "util.h"

using std::auto_ptr;
using std::vector;
using std::endl;
using std::string;
using std::runtime_error;

Output_stream::Output_stream()
{ }

#ifndef _MSC_VER
Output_stream::~Output_stream()
{}
#endif

#ifdef _MSC_VER
Output_stream::Output_stream(const string &file_name) :
	file_name_(file_name),
	f_(file_name.length() == 0 ? stdout : fopen(file_name.c_str(), "wb"))
{
	if (f_ == 0) throw File_open_exception(file_name_);
}
#else
Output_stream::Output_stream(const string &file_name) :
	file_name_(file_name),
	fd_(file_name.length() == 0 ? 1 : open64(file_name.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR))
{
	if (fd_ < 0) {
		perror(0);
		throw File_open_exception(file_name_);
	}
}
#endif

void Output_stream::remove()
{
	if (::remove(file_name_.c_str()) != 0)
		std::cerr << "Warning: Failed to delete file " << file_name_ << std::endl;
}

void Output_stream::close()
{
#ifdef _MSC_VER
	if (f_ && f_ != stdout) {
		fclose(f_);
		f_ = 0;
	}
#else
	if (fd_ != 1 && ::close(fd_) != 0) {
		perror(0);
		throw std::runtime_error(string("Error closing file ") + file_name_);
	}
#endif
}

void Output_stream::write_raw(const char *ptr, size_t count)
{
#ifdef _MSC_VER
	size_t n;
	if ((n = fwrite((const void*)ptr, 1, count, f_)) != count)
		throw File_write_exception(file_name_);
#else
	while (count > 0) {
		const ssize_t n = ::write(fd_, (const void*)ptr, count);
		if (n < 0) {
			perror(0);
			throw std::runtime_error(string("Error writing to file ") + file_name_);
		}
		count -= n;
		ptr += n;
	}
#endif
}

void Output_stream::write_c_str(const string &s)
{
	write(s.c_str(), s.length() + 1);
}

void Output_stream::seekp(size_t p)
{
#ifdef _MSC_VER
	if (_fseeki64(f_, (int64_t)p, SEEK_SET) != 0) throw File_write_exception(file_name_);
#else
	if (lseek64(fd_, p, SEEK_SET) < 0) {
		perror(0);
		throw std::runtime_error("Error calling lseek.");
	}
#endif
}

size_t Output_stream::tell()
{
#ifdef _MSC_VER
	int64_t x;
	if ((x = _ftelli64(f_)) == (int64_t)-1)
		throw std::runtime_error("Error executing ftell on stream " + file_name_);
	return (size_t)x;
#else
	const off64_t n = lseek64(fd_, 0, SEEK_CUR);
	if (n < 0) {
		perror(0);
		throw std::runtime_error("Error calling lseek.");
	}
	return n;
#endif
}

Input_stream::Input_stream(const string &file_name) :
	file_name(file_name),
	line_count(0),
#ifdef _MSC_VER
	f_(file_name.empty() ? stdin : fopen(file_name.c_str(), "rb")),
#else
	fd_(file_name.empty() ? 0 : open64(file_name.c_str(), O_RDONLY)),
#endif
	line_buf_used_(0),
	line_buf_end_(0),
	putback_line_(false),
	eof_(false)
{
#ifdef _MSC_VER
	if (f_ == 0)
		throw File_open_exception(file_name);
#else
	if (fd_ < 0) {
		perror(0);
		throw std::runtime_error(string("Error opening file ") + file_name);
	}
#endif
}

void Input_stream::rewind()
{
#ifdef _MSC_VER
	::rewind(f_);
#else
	seek(0);
#endif
	line_count = 0;
	line_buf_used_ = 0;
	line_buf_end_ = 0;
	putback_line_ = false;
	eof_ = false;
	line.clear();
}

#ifdef _MSC_VER
Input_stream::Input_stream(const Output_stream &tmp_file) :
	file_name(tmp_file.file_name_),
	f_(tmp_file.f_)
{
	::rewind(f_);
}
#else
Input_stream::Input_stream(const Output_stream &tmp_file) :
	file_name(tmp_file.file_name_),
	fd_(tmp_file.fd_)
{
	seek(0);
}
#endif

void Input_stream::seek(size_t pos)
{
#ifdef _MSC_VER
	if (_fseeki64(f_, (int64_t)pos, SEEK_SET) != 0)
		throw std::runtime_error("Error executing seek on file " + file_name);
#else
	if (lseek64(fd_, pos, SEEK_SET) < 0) {
		perror(0);
		throw std::runtime_error("Error calling lseek.");
	}
#endif
}

void Input_stream::seek_forward(size_t n)
{
#ifdef _MSC_VER
	if (_fseeki64(f_, (int64_t)n, SEEK_CUR) != 0)
		throw std::runtime_error("Error executing seek on file " + file_name);
#else
	if (lseek64(fd_, n, SEEK_CUR) < 0) {
		perror(0);
		throw std::runtime_error("Error calling lseek.");
	}
#endif
}

bool Input_stream::eof() const
{
	return eof_;
}

size_t Input_stream::read_bytes(char *ptr, size_t count)
{
#ifdef _MSC_VER
	size_t n;
	if ((n = fread(ptr, 1, count, f_)) != count) {
		if (feof(f_) != 0)
			return n;
		else
			throw File_read_exception(file_name);
	}
	return n;
#else
	ssize_t total = 0;
	while (count > 0) {
		const ssize_t n = ::read(fd_, ptr, count);
		if (n < 0) {
			perror(0);
			throw std::runtime_error(string("Error reading file ") + file_name);
		}
		total += n;
		if (n == 0)
			return total;
		count -= n;
		ptr += n;
	}
	return total;
#endif
}

void Input_stream::read_c_str(string &s)
{
	char c;
	s.clear();
	while (true) {
		if (read(&c, 1) != 1)
			throw std::runtime_error("Unexpected end of file.");
		if (c == 0)
			break;
		s += (char)c;
	}
}

void Input_stream::close_and_delete()
{
	close();
#ifdef _MSC_VER
	if (remove(file_name.c_str()) != 0)
		std::cerr << "Warning: Failed to delete temporary file " << file_name << std::endl;
#endif
}

void Input_stream::close()
{
#ifdef _MSC_VER
	if (f_) {
		fclose(f_);
		f_ = 0;
	}
#else
	if (fd_ != 0 && ::close(fd_) != 0) {
		perror(0);
		throw std::runtime_error(string("Error closing file ") + file_name);
	}
#endif
}

void Input_stream::getline()
{
	if (putback_line_) {
		putback_line_ = false;
		++line_count;
		return;
	}
	line.clear();
	while (true) {
		const char *p = (const char*)memchr(&line_buf_[line_buf_used_], '\n', line_buf_end_ - line_buf_used_);
		if (p == 0) {
			line.append(&line_buf_[line_buf_used_], line_buf_end_ - line_buf_used_);
			line_buf_end_ = read_bytes(line_buf_, line_buf_size);
			line_buf_used_ = 0;
			if (line_buf_end_ == 0) {
				eof_ = true;
				++line_count;
				return;
			}
		}
		else {
			const size_t n = (p - line_buf_) - line_buf_used_;
			line.append(&line_buf_[line_buf_used_], n);
			line_buf_used_ += n + 1;
			const size_t s = line.length() - 1;
			if (!line.empty() && line[s] == '\r')
				line.resize(s);
			++line_count;
			return;
		}
	}
}

void Input_stream::putback_line()
{
	putback_line_ = true;
	--line_count;
}

unsigned Temp_file::n = 0;
uint64_t Temp_file::hash_key;

Temp_file::Temp_file()
{
	if (n == 0) {
#ifdef WIN32
		LARGE_INTEGER count;
		QueryPerformanceCounter(&count);
		hash_key = (uint64_t)(count.HighPart + count.LowPart + count.QuadPart + GetCurrentProcessId());
#else
		timeval count;
		gettimeofday(&count, NULL);
		hash_key = count.tv_sec + count.tv_usec + getpid();
#endif
	}
	std::stringstream ss;
	ss.setf(std::ios::hex, std::ios::basefield);
	if (config.tmpdir != "")
		ss << config.tmpdir << dir_separator;
	ss << "diamond-" << hash_key << "-" << n++ << ".tmp";
	ss >> this->file_name_;
#ifdef _MSC_VER
	this->f_ = fopen(this->file_name_.c_str(), "w+b");
	if (this->f_ == 0)
		throw std::runtime_error("Error opening temporary file: " + this->file_name_);
#else
	this->fd_ = open64(this->file_name_.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
	if (this->fd_ < 0) {
		perror(0);
		throw std::runtime_error(string("Error opening temporary file ") + this->file_name_);
	}
	if (unlink(this->file_name_.c_str()) < 0) {
		perror(0);
		throw std::runtime_error("Error calling unlink.");
	}
#endif
}

string Temp_file::get_temp_dir()
{
	Temp_file t;
	Input_stream f(t);
	f.close_and_delete();
	return extract_dir(f.file_name);
}
