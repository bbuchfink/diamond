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
#include <mutex>
#include <iostream>
#include "../basic/config.h"
#include "log_stream.h"
#include "util.h"
#include "escape_sequences.h"

using namespace std;

Message_stream message_stream;
Message_stream verbose_stream (false);
Message_stream log_stream (false);
std::mutex Message_stream::mtx;

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

Message_stream& Message_stream::operator<<(std::ostream& (*_Pfn)(std::ostream&))
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

Message_stream::Message_stream(bool to_cout, bool to_file) :
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