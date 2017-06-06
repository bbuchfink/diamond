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

#include <vector>
#include <string>
#include <stdexcept>

using std::vector;
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
		runtime_error(string("Error writing file ") + file_name)
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
	FILE *f_;
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

	FILE *f_;
	char line_buf_[line_buf_size];
	size_t line_buf_used_, line_buf_end_;
	bool putback_line_, eof_;

};

#endif /* BINARY_FILE_H_ */
