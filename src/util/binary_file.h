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
	Output_stream();
#ifndef _MSC_VER
	virtual ~Output_stream();
#endif
	Output_stream(const string &file_name);
	void remove();
	virtual void close();
	virtual void write_raw(const char *ptr, size_t count);
	template<typename _t>
	void write(const _t *ptr, size_t count)
	{
		write_raw((const char*)ptr, sizeof(_t)*count);
	}
	template<class _t>
	void write(const vector<_t> &v)
	{
		write(v.data(), v.size());
	}
	void write_c_str(const string &s);
	void seekp(size_t p);
	size_t tell();
protected:
	string file_name_;
#ifdef _MSC_VER
	FILE *f_;
#else
	int fd_;
#endif
	friend struct Input_stream;
};

struct Input_stream
{
	Input_stream(const string &file_name);
	void rewind();
	Input_stream(const Output_stream &tmp_file);
	void seek(size_t pos);
	void seek_forward(size_t n);
	bool eof() const;
	virtual size_t read_bytes(char *ptr, size_t count);
	template<class _t>
	size_t read(_t *ptr, size_t count)
	{
		return read_bytes((char*)ptr, count*sizeof(_t)) / sizeof(_t);
	}
	void read_c_str(string &s);
	void close_and_delete();
	virtual void close();
	void putback(char c);
	void getline();
	void putback_line();

	string file_name;
	string line;
	size_t line_count;
	
protected:

	enum { line_buf_size = 256 };

#ifdef _MSC_VER
	FILE *f_;
#else
	int fd_;
#endif
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
