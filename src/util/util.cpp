/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#include <algorithm>
#include <numeric>
#include <iostream>
#include <mutex>
#include "../basic/config.h"
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

std::vector<std::string> tokenize(const char* str, const char* delimiters)
{
	std::vector<std::string> out;
	std::string token;
	while (*str != 0) {
		while (*str != 0 && strchr(delimiters, *str))
			++str;
		token.clear();
		while (*str != 0 && strchr(delimiters, *str) == nullptr)
			token += *(str++);
		if (token.length() > 0)
			out.push_back(token);
	}
	if (out.size() == 0)
		out.push_back(std::string());
	return out;
}

std::set<int32_t> parse_csv(const std::string& s)
{
	std::set<int32_t> r;
	std::vector<std::string> t(tokenize(s.c_str(), ","));
	for (std::vector<std::string>::const_iterator i = t.begin(); i != t.end(); ++i)
		if (!i->empty()) r.insert(atoi(i->c_str()));
	return r;
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
	if (c < 32 && c >= 0)
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