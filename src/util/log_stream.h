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

#include <iostream>
#include <fstream>
#include "Timer.h"
#include "tinythread.h"

using std::endl;

struct Message_stream
{
	Message_stream(bool to_cout = true, bool to_file = false) :
		to_cout_(to_cout),
		to_file_(to_file)
	{}
	template<typename _t>
	Message_stream& operator<<(const _t& x)
	{
		if(to_cout_)
			std::cout << x;
		if (to_file_) {
			std::ofstream f("diamond.log", std::ios_base::out | std::ios_base::app);
			f << x;
			f.close();
		}
		return *this;
	}
	//Message_stream& operator<<(std::ostream & (__cdecl *_Pfn)(std::ostream&))
	Message_stream& operator<<(std::ostream & (*_Pfn)(std::ostream&))
	{
		if(to_cout_)
			((*_Pfn)(std::cout));
		if (to_file_) {
			mtx.lock();
			std::ofstream f("diamond.log", std::ios_base::out | std::ios_base::app);
			((*_Pfn)(f));
			f.close();
			mtx.unlock();
		}
		return *this;
	}
	static tthread::mutex mtx;
private:
	bool to_cout_, to_file_;
};

extern Message_stream message_stream;
extern Message_stream verbose_stream;
extern Message_stream log_stream;

struct task_timer
{
	task_timer(const char *msg, unsigned level=1) :
		level_(level),
		msg_(msg)
	{
		start(msg);
	}
	~task_timer()
	{
		finish();
	}
	void go(const char *msg)
	{
		finish();
		start(msg);
		msg_ = msg;
	}
	void finish()
	{
		if (!msg_)
			return;
		//if (print_ && !Cfg::debug_log)
		get_stream() << " [" << timer_.getElapsedTimeInSec() << "s]" << endl;
		/*else if (Cfg::debug_log) {
			log_stream << '/' << msg_ << " [" << timer_.getElapsedTimeInSec() << "s]" << endl;
		}*/
		msg_ = 0;
	}
	Message_stream& get_stream() const
	{
		switch (level_) {
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
private:
	void start(const char *msg)
	{
		timer_.start();
		//if (print_ && !Cfg::debug_log) {
			get_stream() << msg << "... " << std::flush;
			//fflush(stdout);
		/**}
		else if (Cfg::debug_log)
			log_stream << msg << "..." << endl;*/
	}
	unsigned level_;
	const char *msg_;
	Timer timer_;
};

#endif
