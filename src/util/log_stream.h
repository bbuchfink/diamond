/****
Copyright � 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <exception>
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
#if defined(__cpp_lib_uncaught_exceptions) && __cpp_lib_uncaught_exceptions >= 201411L
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
