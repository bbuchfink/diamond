/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <exception>
#include <ostream>
#include <chrono>
#include <stdint.h>
#include <limits.h>

extern std::ostream* message_stream;
extern std::ostream* log_stream;

void cleanup();

struct TaskTimer
{
	TaskTimer(std::ostream* stream, unsigned level = 1) :
		level_(level),
		msg_(nullptr),
		stream_(stream)
	{
		start(nullptr);
	}
	TaskTimer(unsigned level = 1) :
		TaskTimer(get_stream(level), level)
	{}
	TaskTimer(const char* msg, std::ostream* stream, unsigned level = 1) :
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
		*stream_ << " [" << get() << "s]" << std::endl;
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
		*stream_ << msg << "... " << std::flush;

	}
	std::ostream* get_stream(unsigned level) const
	{
		switch (level) {
		case 1:
			return message_stream;
		case 2:
		case 3:
			return log_stream;
		default:
			return message_stream;
		}
	}
	unsigned level_;
	const char *msg_;
	std::ostream* stream_;
	std::chrono::high_resolution_clock::time_point t;
};