/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef BINARY_FILE_H_
#define BINARY_FILE_H_

#include <memory>
#include <vector>
#include <string>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <assert.h>
#include "log_stream.h"
#include "system.h"

using std::auto_ptr;
using std::vector;
using std::endl;
using std::string;
using std::runtime_error;

struct File_read_exception : public std::runtime_error
{
	File_read_exception(const string &file_name) :
		runtime_error(string("Error reading file ") + file_name)
	{ }
};

struct File_write_exception : public std::runtime_error
{
	File_write_exception(const string &file_name) :
		runtime_error(string("Error writing file ") + file_name + ". Disk full?")
	{ }
};

struct File_open_exception : public std::runtime_error
{
	File_open_exception(const string &file_name) :
		std::runtime_error(string("Error opening file " + file_name))
	{ }
};

struct Output_stream
{

	Output_stream()
	{ }
#ifndef _MSC_VER
	virtual ~Output_stream()
	{}
#endif
	Output_stream(const string &file_name) :
		file_name_(file_name),
		f_(file_name.length() == 0 ? stdout : fopen(file_name.c_str(), "wb"))
	{
		if (f_ == 0) throw File_open_exception(file_name_);
	}
	void remove()
	{
		if (::remove(file_name_.c_str()) != 0)
			std::cerr << "Warning: Failed to delete file " << file_name_ << std::endl;
	}
	virtual void close()
	{
		if (f_ && f_ != stdout) {
			fclose(f_);
			f_ = 0;
		}
	}
	virtual void write(const char *ptr, size_t count)
	{
		size_t n;
		if ((n = fwrite((const void*)ptr, 1, count, f_)) != count)
			throw File_write_exception(file_name_);
	}
	template<typename _t>
	void typed_write(const _t *ptr, size_t count)
	{
		size_t n;
		if ((n = fwrite((const void*)ptr, sizeof(_t), count, f_)) != count)
			throw File_write_exception(file_name_);
	}
	template<class _t>
	void write(const vector<_t> &v, bool write_size = true)
	{
		size_t size = v.size();
		if (write_size)
			typed_write(&size, 1);
		typed_write(v.data(), size);
	}
	void write_c_str(const string &s)
	{
		write(s.c_str(), s.length() + 1);
	}
	void seekp(size_t p)
	{
		if (FSEEK(f_, (int64_t)p, SEEK_SET) != 0) throw File_write_exception(file_name_);
	}
	size_t tell()
	{
		int64_t x;
		if ((x = FTELL(f_)) == (int64_t)-1)
			throw std::runtime_error("Error executing ftell on stream " + file_name_);
		return (size_t)x;
	}
protected:
	string file_name_;
	FILE *f_;
	friend struct Input_stream;
};

struct Input_stream
{

	Input_stream(const string &file_name):
		file_name (file_name),
		line_count(0),
		f_(fopen(file_name.c_str(), "rb")),		
		line_buf_used_(0),
		line_buf_end_(0),
		putback_line_(false),
		eof_(false)
	{
		if (f_ == 0)
			throw File_open_exception(file_name);
	}

	void rewind()
	{
		::rewind(f_);
		line_count = 0;
		line_buf_used_ = 0;
		line_buf_end_ = 0;
		putback_line_ = false;
		eof_ = false;
		line.clear();
	}

	Input_stream(const Output_stream &tmp_file):
		file_name (tmp_file.file_name_),
		f_(tmp_file.f_)
	{
		::rewind(f_);
	}

	void seek(size_t pos)
	{
		if (FSEEK(f_, (int64_t)pos, SEEK_SET) != 0)
			throw std::runtime_error("Error executing seek on file " + file_name);
	}

	void seek_forward(size_t n)
	{
		if (FSEEK(f_, (int64_t)n, SEEK_CUR) != 0)
			throw std::runtime_error("Error executing seek on file " + file_name);
	}

	bool eof() const
	{
		return eof_;
	}

	virtual size_t read_bytes(char *ptr, size_t count)
	{
		size_t n;
		if ((n = fread(ptr, 1, count, f_)) != count) {
			if (feof(f_) != 0)
				return n;
			else
				throw File_read_exception(file_name);
		}
		return n;
	}

	template<class _t>
	size_t read(_t *ptr, size_t count)
	{
		return read_bytes((char*)ptr, count*sizeof(_t)) / sizeof(_t);
	}

	template<class _t>
	void read(vector<_t> &v)
	{
		size_t size;
		if (read(&size, 1) != 1)
			throw File_read_exception(file_name);
		v.clear();
		v.resize(size);
		if (read(v.data(), size) != size)
			throw File_read_exception(file_name);
	}

	template<class _t>
	void skip_vector()
	{
		size_t size;
		if (read(&size, 1) != 1)
			throw File_read_exception(file_name);
		seek_forward(size * sizeof(_t));
	}

	void read_c_str(string &s)
	{
		int c;
		s.clear();
		while ((c = fgetc(f_)) != 0) {
			if (c == EOF)
				throw File_read_exception(file_name);
			s += (char)c;
		}
	}

	void close_and_delete()
	{
		close();
#ifdef _MSC_VER
		if (remove(file_name.c_str()) != 0)
			std::cerr << "Warning: Failed to delete temporary file " << file_name << std::endl;
#endif
	}

	virtual void close()
	{
		if (f_) {
			fclose(f_);
			f_ = 0;
		}
	}

	void putback(char c)
	{
		if (ungetc(c, f_) != (int)c)
			throw File_read_exception(file_name);
	}

	void getline()
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

	void putback_line()
	{
		putback_line_ = true;
		--line_count;
	}

	string file_name;
	string line;
	size_t line_count;
	
protected:

	enum { line_buf_size = 256 };

	FILE *f_;
	char line_buf_[line_buf_size];
	size_t line_buf_used_, line_buf_end_;
	bool putback_line_, eof_;

};

struct Buffered_file : public Input_stream
{

	Buffered_file(const string& file_name):
		Input_stream (file_name)
	{ init(); }

	Buffered_file(const Output_stream &tmp_file):
		Input_stream (tmp_file)
	{ init(); }

	bool eof()
	{
		if(ptr_ < end_)
			return false;
		else if(end_ < &data_[buffer_size])
			return true;
		else {
			fetch();
			return eof();
		}
	}

	template<typename _t>
	void read(_t &dst)
	{
		const char *const p = ptr_ + sizeof(_t);
		if(p > end_) {
			if(end_ < &data_[buffer_size])
				throw File_read_exception(file_name);
			fetch();
			return read(dst);
		}
		dst = *reinterpret_cast<_t*>(ptr_);
		ptr_ += sizeof(_t);
	}

	template<typename _t>
	void read(_t* dst, size_t n)
	{
		for(size_t i=0;i<n;++i)
			read(*(dst++));
	}

	void read_c_str(string &dst)
	{
		dst.clear();
		char c;
		while(read(c), c != '\0')
			dst.push_back(c);
	}

	void read_packed(uint8_t length, uint32_t &dst)
	{
		switch(length) {
		case 0: uint8_t x; read(x); dst = x; break;
		case 1: uint16_t y; read(y); dst = y; break;
		case 2: read(dst);
		}
	}

	const char* ptr() const
	{ return ptr_; }

	void seek(size_t pos)
	{
		//this->clear();
		Input_stream::seek(pos);
		init();
	}

private:

	void fetch()
	{
		ptrdiff_t d = end_ - ptr_;
		assert(d >= 0);
		memcpy(&data_[0], ptr_, d);
		ptr_ = &data_[0];
		end_ = read_block(ptr_+d, buffer_size - d);
	}

	void init()
	{
		ptr_ = &data_[0];
		end_ = read_block(ptr_, buffer_size);
	}

	char* read_block(char* ptr, size_t size)
	{ return ptr + Input_stream::read(ptr, size); }

	enum { buffer_size = 8192 };

	char data_[buffer_size];
	char *ptr_, *end_;

};

struct Buffered_ostream
{

	Buffered_ostream(Output_stream &s):
		stream_ (s),
		ptr_ (&data_[0]),
		end_ (&data_[buffer_size])
	{ }

	template<typename _t>
	void write(const _t &x)
	{
		const char *const p = ptr_ + sizeof(_t);
		if(p > end_) {
			flush();
			return write(x);
		}
		*reinterpret_cast<_t*>(ptr_) = x;
		ptr_ += sizeof(_t);
	}

	void close()
	{
		flush();
		((Output_stream*)this)->close(); 
	}

private:

	void flush()
	{
		stream_.write(data_, ptr_ - data_);
		ptr_ = data_;
	}

	enum { buffer_size = 4096 };

	Output_stream &stream_;
	char data_[buffer_size];
	char *ptr_, * const end_;

};

#endif /* BINARY_FILE_H_ */
