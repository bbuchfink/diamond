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

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/thread/mutex.hpp>

using std::exception;
using std::string;

#define THROW_EXCEPTION(exception, p1) throw exception(p1, __PRETTY_FUNCTION__, __LINE__)

struct diamond_exception : public exception
{
	string msg;
	diamond_exception(const char *function, unsigned line):
		msg(string("function ") + function + " line " + boost::lexical_cast<string>(line))
	{ }
	diamond_exception(const string &msg):
		msg (msg)
	{ }
	~diamond_exception() throw() { }
	virtual const char* what() const throw()
	{ return msg.c_str(); }
};

struct file_io_exception : public diamond_exception
{
	file_io_exception(const string &file_name, const char* function, unsigned line):
		diamond_exception(function, line)
	{ msg += ". Error reading file " + file_name; }
};

struct file_open_exception: public diamond_exception
{
	file_open_exception(const string &file_name, const char* function, unsigned line):
		diamond_exception(function, line)
	{ msg += ". Error opening file " + file_name; }
};

struct file_io_write_exception: public diamond_exception
{
	file_io_write_exception(const string &file_name, const char* function, unsigned line):
		diamond_exception(function, line)
	{ msg += ". Error writing file " + file_name; }
};

struct memory_alloc_exception: public exception
{
	virtual const char* what() const throw()
	{ return "Failed to allocate memory"; }
};

struct hash_table_overflow_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Hash table overflow"; }
};

struct invalid_database_version_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Incompatible database version"; }
};

struct invalid_parameter_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Invalid parameter"; }
};

struct invalid_sequence_char_exception : public exception
{
	const string msg;
	invalid_sequence_char_exception(char ch):
		msg (string("Invalid character (")+ch+"/" + boost::lexical_cast<string>(int(ch)) + ") in sequence")
	{ }
	~invalid_sequence_char_exception() throw()
	{ }
	virtual const char* what() const throw()
	{ return msg.c_str(); }
};

struct file_format_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Invalid input file format"; }
};

struct Exception_state
{
	Exception_state():
		active_ (false)
	{ }
	void set(const std::exception &e)
	{
		if(!active_) {
			mtx_.lock();
			what_ = string(typeid(e).name()) + ": " + e.what();
			active_ = true;
			mtx_.unlock();
		}
	}
	void sync() const
	{
		if(active_)
			throw std::runtime_error(what_);
	}
	bool operator()() const
	{ return active_; }
private:
	boost::mutex mtx_;
	bool active_;
	string what_;
} exception_state;

#endif /* EXCEPTIONS_H_ */
