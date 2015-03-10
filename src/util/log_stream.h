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

#ifndef LOG_STREAM_H_
#define LOG_STREAM_H_

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/timer/timer.hpp>

using std::string;
using std::endl;

boost::iostreams::filtering_ostream verbose_stream;
boost::iostreams::filtering_ostream log_stream;

struct task_timer : public boost::timer::cpu_timer
{
	task_timer(const char *msg, bool print):
		print_ (print),
		msg_ (msg)
	{ start(msg); }
	~task_timer()
	{ finish(); }
	void go(const char *msg)
	{
		finish();
		boost::timer::cpu_timer::start();
		start(msg);
		msg_ = msg;
	}
	void finish()
	{
		if(!msg_)
			return;
		if(print_ && !program_options::debug_log)
			verbose_stream << boost::timer::format(elapsed(), 1, "[%ws]") << endl;
		else {
			log_stream << '/' << msg_ << boost::timer::format(elapsed(), 1, " [%ws]") << endl;
		}
		msg_ = 0;
	}
private:
	void start(const char *msg)
	{
		if(print_ && !program_options::debug_log) {
			verbose_stream << msg << "... " << std::flush;
			fflush(stdout);
		} else
			log_stream << msg << "..." << endl;
	}
	bool print_;
	const char *msg_;
};

#endif
