/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
						Benjamin Buchfink
						Eberhard Karls Universitaet Tuebingen

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

#pragma once
#include <ostream>
#include <fstream>
#include <mutex>
#include <chrono>
#include <stdint.h>
#include <limits.h>

struct MessageStream
{
	MessageStream(bool to_cout = true, bool to_file = false);
	template<typename T>
	MessageStream& operator<<(const T& x)
	{
		if(to_cout_)
			(*out_stream_) << x;
		if (to_file_) {
			std::ofstream f("diamond.log", std::ios_base::out | std::ios_base::app);
			f << x;
			f.close();
		}
		return *this;
	}
	//Message_stream& operator<<(std::ostream & (__cdecl *_Pfn)(std::ostream&))
	MessageStream& operator<<(std::ostream& (*_Pfn)(std::ostream&));
	static std::mutex mtx;
private:
	std::ostream* out_stream_;
	bool to_cout_, to_file_;
};

extern MessageStream message_stream;
extern MessageStream verbose_stream;
extern MessageStream log_stream;

struct TaskTimer
{
	TaskTimer(MessageStream& stream, unsigned level = 1) :
		level_(level),
		msg_(nullptr),
		stream_(stream)
	{
		start(nullptr);
	}
	TaskTimer(unsigned level = 1) :
		TaskTimer(get_stream(level), level)
	{}
	TaskTimer(const char* msg, MessageStream& stream, unsigned level = 1) :
		level_(level),
		msg_(msg),
		stream_(stream)
	{
		start(msg);
	}
	TaskTimer(const char* msg, unsigned level = 1):
		TaskTimer(msg, get_stream(level), level)
	{}
	~TaskTimer()
	{
#if (__cplusplus >= 201703L && !defined(__APPLE__)) || defined(_MSC_VER)
		if (!std::uncaught_exceptions())
#else
		if (!std::uncaught_exception())
#endif
			finish();
	}
	void go(const char* msg = nullptr)
	{
		finish();
		start(msg);
		msg_ = msg;
	}
	void go(const std::string& s) {
		go(s.c_str());
	}
	void finish()
	{
		if (!msg_ || level_ == UINT_MAX)
			return;
		stream_ << " [" << get() << "s]" << std::endl;
		msg_ = 0;
	}
	double get()
	{
		return (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count() / 1000.0;
	}
	int64_t seconds() const {
		return (int64_t)std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t).count();
	}
	int64_t milliseconds() const {
		return (int64_t)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count();
	}
	int64_t microseconds() const {
		return (int64_t)std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t).count();
	}
	int64_t nanoseconds() const {
		return (int64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t).count();
	}
private:
	void start(const char *msg)
	{
		t = std::chrono::high_resolution_clock::now();
		if (level_ == UINT_MAX)
			return;
		if (!msg)
			return;
		stream_ << msg << "... " << std::flush;

	}
	MessageStream& get_stream(unsigned level) const
	{
		switch (level) {
		case 1:
			return message_stream;
		case 2:
			return verbose_stream;
		case 3:
			return log_stream;
		default:
			return message_stream;
		}
	}
	unsigned level_;
	const char *msg_;
	MessageStream& stream_;
	std::chrono::high_resolution_clock::time_point t;
};

void exit_with_error(const std::exception& e);
