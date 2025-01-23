#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <limits.h>
#include <stdint.h>
#include <limits>
#include "string.h"

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

}}
