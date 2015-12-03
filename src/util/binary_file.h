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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "../basic/exceptions.h"
#include "temp_file.h"

namespace io = boost::iostreams;

using std::auto_ptr;
using std::vector;
using std::endl;

struct File_read_exception : public diamond_exception
{
	File_read_exception(const char* file_name, size_t count, ssize_t n):
		diamond_exception (string("Error reading file ") + file_name + " (" + boost::lexical_cast<string>(count) + '/' + boost::lexical_cast<string>(n) + ')')
	{ }
};

struct File_write_exception : public diamond_exception
{
	File_write_exception(const char* file_name, size_t count, ssize_t n):
		diamond_exception (string("Error writing file ") + file_name + " (" + boost::lexical_cast<string>(count) + '/' + boost::lexical_cast<string>(n) + ')')
	{ }
};

struct Input_stream : public io::filtering_istream
{

	Input_stream(const string &file_name, bool gzipped = false):
		file_name (file_name)
	{
		if(gzipped) {
			FILE *f = fopen(file_name.c_str(), "rb");
			if(f == 0)
				THROW_EXCEPTION(file_open_exception, file_name);
			unsigned char id[2];
			if(fread(id, 1, 2, f) != 2)
				THROW_EXCEPTION(file_io_exception, file_name);
			fclose(f);
			if(id[0] == 0x1f && id[1] == 0x8b) {
				log_stream << "Detected gzip compressed file " << file_name << endl;
				this->push(io::gzip_decompressor ());
			}
		}
		io::file_source f (file_name, std::ios_base::in | std::ios_base::binary);
		if(!f.is_open())
			THROW_EXCEPTION(file_open_exception, file_name);
		this->push(f);
	}

	Input_stream(const Temp_file &tmp_file)
	{
		if(lseek(tmp_file.file_descriptor(), 0, SEEK_SET) == -1)
			throw std::runtime_error("Error executing seek on temporary file.");
		this->push(io::file_descriptor_source(tmp_file.file_descriptor(), io::close_handle));
	}

	void seek(size_t pos)
	{
		this->seekg(pos);
		if(!this->good())
			THROW_EXCEPTION(file_io_exception, file_name);
	}

	void seek_forward(size_t n)
	{
		this->seekg(n, ios_base::cur);
		if(!this->good())
			THROW_EXCEPTION(file_io_exception, file_name);
	}

	template<class _t>
	size_t read(_t *ptr, size_t count)
	{
		ssize_t n;
		if((n = io::read(*this, reinterpret_cast<char*>(ptr), sizeof(_t) * count)) != (ssize_t)(count * sizeof(_t))) {
			if(n == EOF)
				return 0;
			else if(n >= 0 && this->get() == EOF && (n % sizeof(_t) == 0))
				return n/sizeof(_t);
			else
				throw File_read_exception(file_name.c_str(), sizeof(_t) * count, n);
		}
		return n/sizeof(_t);
	}

	template<class _t>
	void read(vector<_t> &v)
	{
		size_t size;
		if(read(&size, 1) != 1)
			THROW_EXCEPTION(file_io_exception, file_name);
		v.resize(size);
		if(read(v.data(), size) != size)
			THROW_EXCEPTION(file_io_exception, file_name);
	}

	template<class _t>
	void skip_vector()
	{
		size_t size;
		if(read(&size, 1) != 1)
			THROW_EXCEPTION(file_io_exception, file_name);
		seek_forward(size * sizeof(_t));
	}

	void read_c_str(string &s)
	{
		std::getline(*this, s, '\0');
		if(!this->good())
			THROW_EXCEPTION(file_io_exception, file_name);
	}

	void close()
	{
		this->set_auto_close(true);
		this->pop();
	}

	const string file_name;

};

struct Output_stream2
{
	Output_stream2()
	{ }
	Output_stream2(const string &file_name):
		file_name_ (file_name),
		f_ (fopen(file_name.c_str(), "wb"))
	{ if(f_ == 0) THROW_EXCEPTION(file_open_exception, file_name_); }
	void close()
	{ if(f_) fclose(f_); f_ = 0; }
	~Output_stream2()
	{ close(); }
	template<typename _t>
	void write(const _t *ptr, size_t count)
	{
		size_t n;
		if((n=fwrite((const void*)ptr, sizeof(_t), count, f_)) != count)
			throw File_write_exception (file_name_.c_str(), count, n);
	}
	template<class _t>
	void write(const vector<_t> &v)
	{
		size_t size = v.size();
		write(&size, 1);
		write(v.data(), size);
	}
	void seekp(size_t p)
	{ if(fseek(f_, p, SEEK_SET) != 0) throw File_write_exception (file_name_.c_str(), 0, 0); }
private:
	const string file_name_;
	FILE *f_;
};

struct Output_stream : public io::filtering_ostream
{

	enum Flags { file_sink, stdout_sink };

	Output_stream()
	{ }

	Output_stream(const string &file_name, bool gzipped = false, Flags flags = file_sink):
		file_name_ (file_name)
	{
		if(gzipped)
			this->push(io::gzip_compressor ());
		if(flags == file_sink) {
			io::file_sink f (file_name, std::ios_base::out | std::ios_base::binary);
			if(!f.is_open())
				THROW_EXCEPTION(file_open_exception, file_name_);
			this->push(f);
		} else
			this->push(std::cout);
	}

	Output_stream(const Temp_file &tmp_file)
	{
		this->push(io::file_descriptor_sink(tmp_file.file_descriptor(), io::never_close_handle));
	}

	template<typename _t>
	void write(const _t *ptr, size_t count)
	{
		ssize_t n;
		if((n = io::write(*this, reinterpret_cast<const char*>(ptr), sizeof(_t) * count)) != (ssize_t)(sizeof(_t) * count))
			throw File_write_exception (file_name_.c_str(), sizeof(_t) * count, n);
	}

	void write_c_str(const string &s)
	{ write(s.c_str(), s.length()+1); }

	template<class _t>
	void write(const vector<_t> &v, bool write_size = true)
	{
		size_t size = v.size();
		if(write_size)
			write(&size, 1);
		write(v.data(), size);
	}

	size_t tell()
	{
		if(this->tellp() < 0)
			throw std::runtime_error("Unable to execute tellp() on output stream.");
		return this->tellp();
	}

	void seek(size_t pos)
	{
		this->seekp(pos);
		if(!this->good())
			throw std::runtime_error("Unable to execute seekp() on output stream.");
	}

	void close()
	{
		this->set_auto_close(true);
		this->pop();
	}

private:

	const string file_name_;

};

struct Buffered_file : public Input_stream
{

	Buffered_file(const string& file_name, bool gzipped=false):
		Input_stream (file_name, gzipped)
	{ init(); }

	Buffered_file(const Temp_file &tmp_file):
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
				THROW_EXCEPTION(file_io_exception, file_name);
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
		this->clear();
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

	~Buffered_ostream()
	{ flush(); }

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

void copy_file(boost::iostreams::filtering_ostream &s, const char* file_name)
{
	s << file_name << endl;
	try {
		Input_stream f (file_name);
		s << f.rdbuf() << endl;
	} catch (file_open_exception &e) {
	}
}

#endif /* BINARY_FILE_H_ */
