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

#include <numeric>
#include <iostream>
#include "log_stream.h"
#include "util.h"
#include "escape_sequences.h"

using std::string;
using std::vector;
using std::endl;

class Nullbuf : public std::streambuf {
protected:
	int overflow(int c) override {
		return c;
	}
};

class Nullostream : public std::ostream {
public:
	Nullostream() : std::ostream(&buf_) {}
private:
	Nullbuf buf_;
};

std::ostream* message_stream = new Nullostream;
std::ostream* log_stream = new Nullostream;

void cleanup() {
	if (message_stream != &std::cerr)
		delete message_stream;
	if (log_stream != &std::cerr)
		delete log_stream;
}

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

void print_binary(uint64_t x)
{
	for (unsigned i = 0; i < 64; ++i) {
		std::cout << (x & 1);
		x >>= 1;
	}
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