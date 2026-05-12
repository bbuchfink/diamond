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

#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <limits.h>
#include <stdint.h>
#include <limits>
#include "string.h"
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <algorithm>


using std::string;
using std::to_string;
using std::stringstream;
using std::setprecision;
using std::fixed;
using std::runtime_error;
using std::set;
using std::vector;
using std::numeric_limits;

string convert_size(size_t size) {
	static const char *SIZES[] = { "B", "KB", "MB", "GB", "TB", "PB" };
	size_t div = 0;
	size_t rem = 0;

	while (size >= 1024 && div < (sizeof SIZES / sizeof *SIZES)) {
		rem = (size % 1024);
		div++;
		size /= 1024;
	}

	stringstream ss;
	ss << std::fixed << std::setprecision(1) << (double)size + (double)rem / 1024.0 << ' ' << SIZES[div];
	return ss.str();
}

namespace Util { namespace String {

string replace(const string& s, char a, char b) {
	string r = s;
	size_t i = r.find_first_of(a);
	while (i != string::npos) {
		r[i++] = b;
		i = r.find_first_of(a, i);
	}
	return r;
}

std::string ratio_percentage(const double x, const double y) {
	stringstream ss;
	ss << fixed << setprecision(0) << x << '/' << y << " (" << setprecision(2) << x / y * 100.0 << "%)";
	return ss.str();
}

std::string ratio_percentage(const size_t x, const size_t y) {
	return ratio_percentage((double)x, (double)y);
}

int64_t interpret_number(const std::string& s) {
	stringstream ss(s);
	double n;
	ss >> n;
	char c = 0;
	ss >> c;
	if (ss.eof())
		throw runtime_error("Missing size specifier in number: " + s + ". Permitted values: T, G, M, K");
	double mult = 1;
	switch(c) {
		case 'T':
		case 't':
			mult = 1e12;
			break;
		case 'G':
		case 'g':
			mult = 1e9;
			break;
		case 'M':
		case 'm':
			mult = 1e6;
			break;
		case 'K':
		case 'k':
			mult = 1e3;
			break;
		default:
			throw runtime_error(string("Invalid size specifier (") + c + ") in number: " + s + ". Permitted suffixes: T, G, M, K");
			break;
	}
	ss >> c;
	if (!ss.eof())
		throw runtime_error("Invalid number format: " + s);
	return int64_t(n * mult);
}

vector<string> tokenize(const char* str, const char* delimiters)
{
	vector<string> out;
	string token;
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
		out.push_back(string());
	return out;
}

set<int32_t> parse_csv(const string& s)
{
	set<int32_t> r;
	vector<string> t(tokenize(s.c_str(), ","));
	for (vector<string>::const_iterator i = t.begin(); i != t.end(); ++i)
		if (!i->empty()) r.insert(atoi(i->c_str()));
	return r;
}

template<>
int64_t convert_string<int64_t>(const char* s) {
	char* end;
	long long i = strtoll(s, &end, 10);
	if ((i == 0 && strcmp(s, "0") != 0) || i == LLONG_MAX || i == LLONG_MIN || *end != '\0' || i < numeric_limits<int64_t>::min() || i > numeric_limits<int64_t>::max())
		throw runtime_error(string("Error converting integer value: ") + s);
	return i;
}

template<>
int32_t convert_string<int32_t>(const char* s) {
	const int64_t i = convert_string<int64_t>(s);
	if (i < (int64_t)INT32_MIN || i >(int64_t)INT32_MAX)
		throw runtime_error(string("Error converting integer value: ") + s);
	return (int32_t)i;
}

template<>
uint64_t convert_string<uint64_t>(const char* s) {
	char* end;
	unsigned long long i = strtoull(s, &end, 10);
	if ((i == 0 && strcmp(s, "0") != 0) || i == ULLONG_MAX || *end != '\0' || i > numeric_limits<uint64_t>::max())
		throw runtime_error(string("Error converting integer value: ") + s);
	return i;
}

template<>
uint32_t convert_string<uint32_t>(const char* s) {
	const int64_t i = convert_string<int64_t>(s);
	if (i < 0 || i > (int64_t)UINT32_MAX)
		throw runtime_error(string("Error converting integer value: ") + s);
	return (uint32_t)i;
}

string format(double number) {
	const vector<string> suffixes = { "", "K", "M", "G", "T", "P", "E" };

	if (number == 0) {
		return "0";
	}

	bool is_negative = number < 0;
	double abs_number = std::abs(number);

	int suffix_index = 0;
	if (abs_number >= 1000) {
		int exp = static_cast<int>(std::log10(abs_number) / 3);
		suffix_index = std::min(exp, static_cast<int>(suffixes.size() - 1));
	}

	double divisor = std::pow(1000.0, suffix_index);
	double formatted_number = abs_number / divisor;

	std::ostringstream oss;
	oss << std::fixed << std::setprecision(2) << formatted_number;
	string num_str = oss.str();

	// Trim trailing zeros and decimal point if needed
	size_t dot_pos = num_str.find('.');
	if (dot_pos != std::string::npos) {
		// Find the last non-zero character after the decimal point
		size_t last_non_zero = num_str.find_last_not_of('0');
		if (last_non_zero != std::string::npos) {
			num_str = num_str.substr(0, last_non_zero + 1);
			// If the decimal point is now at the end, remove it
			if (num_str.back() == '.') {
				num_str.pop_back();
			}
		}
	}

	// Check if after trimming, the number is 1000 or more, adjust suffix if possible
	if (formatted_number >= 999.995 && suffix_index < (int)suffixes.size() - 1) {
		suffix_index++;
		divisor *= 1000;
		formatted_number = abs_number / divisor;
		oss.str("");
		oss << std::fixed << std::setprecision(2) << formatted_number;
		num_str = oss.str();
		// Trim again
		dot_pos = num_str.find('.');
		if (dot_pos != std::string::npos) {
			size_t last_non_zero = num_str.find_last_not_of('0');
			if (last_non_zero != std::string::npos) {
				num_str = num_str.substr(0, last_non_zero + 1);
				if (num_str.back() == '.') {
					num_str.pop_back();
				}
			}
		}
	}

	string result = (is_negative ? "-" : "") + num_str + suffixes[suffix_index];
	return result;
}

string format(unsigned long number) {
	return format((double)number);
}

string format(unsigned long long number) {
	return format((double)number);
}

string format(long number) {
	return format((double)number);
}

string format(long long number) {
	return format((double)number);
}

}}