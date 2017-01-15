/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include <stdio.h>
#include <exception>
#include <sstream>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include "../basic/config.h"
#include "log_stream.h"
#include "complexity_filter.h"
#include "util.h"
#include "temp_file.h"
#include "binary_file.h"

Message_stream message_stream;
Message_stream verbose_stream (false);
Message_stream log_stream (false);
tthread::mutex Message_stream::mtx;

const Complexity_filter Complexity_filter::instance;
TLS_PTR vector<Ptr_wrapper_base*> *TLS::ptr_;

#ifndef _MSC_VER
const char dir_separator = '/';
#else
const char dir_separator = '\\';
#endif

string extract_dir(const string & s)
{
	return s.find_last_of(dir_separator) == string::npos ? "" : s.substr(0, s.find_last_of(dir_separator));
}

unsigned Temp_file::n = 0;
uint64_t Temp_file::hash_key;

Temp_file::Temp_file()
{
	if (n == 0) {
#ifdef WIN32
		LARGE_INTEGER count;
		QueryPerformanceCounter(&count);
		hash_key = (uint64_t)(count.HighPart + count.LowPart + count.QuadPart + GetCurrentProcessId());
#else
		timeval count;
		gettimeofday(&count, NULL);
		hash_key = count.tv_sec + count.tv_usec + getpid();
#endif
	}
	std::stringstream ss;
	ss.setf(std::ios::hex, std::ios::basefield);
	if (config.tmpdir != "")
		ss << config.tmpdir << dir_separator;
	ss << "diamond-" << hash_key << "-" << n++ << ".tmp";
	ss >> this->file_name_;
	this->f_ = fopen(this->file_name_.c_str(), "w+b");
	if (this->f_ == 0)
		throw std::runtime_error("Error opening temporary file: " + this->file_name_);
#ifndef _MSC_VER
	unlink(this->file_name_.c_str());
#endif
}

string Temp_file::get_temp_dir()
{
	Temp_file t;
	Input_stream f(t);
	f.close_and_delete();
	return extract_dir(f.file_name);
}

Sd::Sd(const vector<Sd> &groups):
	A(0),
	Q(0),
	k(0)
{
	for (unsigned i = 0; i < groups.size(); ++i) {
		k += groups[i].k;
		A += groups[i].A * groups[i].k;
		Q += groups[i].Q;
	}
	A /= k;
	for (unsigned i = 0; i < groups.size(); ++i)
		Q += (groups[i].A - A)*(groups[i].A - A)*groups[i].k;
}