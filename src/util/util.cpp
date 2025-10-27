/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

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

#include <numeric>
#include <iostream>
#include <mutex>
#include "log_stream.h"
#include "util.h"
#include "escape_sequences.h"
#include "profiler.h"

using std::string;
using std::vector;
using std::endl;

MessageStream message_stream;
MessageStream verbose_stream (false);
MessageStream log_stream (false);
std::mutex MessageStream::mtx;
std::map<std::string, uint64_t> Profiler::times;

#ifndef _MSC_VER
const char dir_separator = '/';
#else
const char dir_separator = '\\';
#endif

string extract_dir(const string & s)
{
	return s.find_last_of(dir_separator) == string::npos ? "" : s.substr(0, s.find_last_of(dir_separator));
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

const EscapeSequence EscapeSequences::xml_data[5] = {
	{ '\"', "&quot;" },
	{ '\'', "&apos;" },
	{ '<', "&lt;" },
	{ '>', "&gt;" },
	{ '&', "&amp;" }
};

const EscapeSequences EscapeSequences::XML(EscapeSequences::xml_data, 5);

string join(const char *c, const vector<string> &v) {
	string s;
	if (v.empty())
		return s;
	s.reserve(accumulate(v.begin(), v.end(), (size_t)0, [](size_t a, const string &b) { return a + b.length(); }) + v.size() - 1);
	for (size_t i = 0; i < v.size() - 1; ++i) {
		s += v[i];
		s += c;
	}
	s += v.back();
	return s;
}

MessageStream& MessageStream::operator<<(std::ostream& (*_Pfn)(std::ostream&))
{
	if (to_cout_)
		((*_Pfn)(std::cerr));
	if (to_file_) {
		mtx.lock();
		std::ofstream f("diamond.log", std::ios_base::out | std::ios_base::app);
		((*_Pfn)(f));
		f.close();
		mtx.unlock();
	}
	return *this;
}

MessageStream::MessageStream(bool to_cout, bool to_file) :
	out_stream_(&std::cerr),
	to_cout_(to_cout),
	to_file_(to_file)	
{}

void print_binary(uint64_t x)
{
	for (unsigned i = 0; i < 64; ++i) {
		std::cout << (x & 1);
		x >>= 1;
	}
}

void exit_with_error(const std::exception& e) {
	std::cerr << "Error: " << e.what() << endl;
	log_stream << "Error: " << e.what() << endl;
	exit(EXIT_FAILURE);
}

std::string to_upper_case(const std::string& s)
{
	std::string r;
	for (std::string::const_iterator i = s.begin(); i != s.end(); ++i)
		r.push_back(toupper(*i));
	return r;
}

std::string to_lower_case(const std::string& s)
{
	std::string r;
	for (std::string::const_iterator i = s.begin(); i != s.end(); ++i)
		r.push_back(tolower(*i));
	return r;
}

std::string print_char(char c)
{
	char buf[16];
	if (c < 32)
		snprintf(buf, sizeof(buf), "ASCII %u", (unsigned)c);
	else
		snprintf(buf, sizeof(buf), "%c", c);
	return std::string(buf);
}

std::string hex_print(const char* x, int len) {
	std::string out;
	char d[3];
	for (int i = 0; i < len; i++) {
		snprintf(d, sizeof(d), "%02x", (unsigned char)x[i]);
		out += d;
	}
	return out;
}