#pragma once
#include <string>
#include <string.h>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <cstdint>
#include <set>

inline bool ends_with(const std::string &s, const char *t) {
	if (s.length() < strlen(t))
		return false;
	return s.compare(s.length() - strlen(t), std::string::npos, t) == 0;
}

inline std::string rstrip(const std::string &s, const char *t) {
	const size_t l = strlen(t);
	if (s.length() < l)
		return s;
	if (s.compare(s.length() - l, std::string::npos, t) == 0)
		return std::string(s.begin(), s.end() - l);
	else
		return s;
}

inline size_t max_len(const char** s, size_t n) {
	size_t l = 0;
	for (size_t i = 0; i < n; ++i)
		l = std::max(l, strlen(s[i]));
	return l;
}

template<typename _it>
inline std::vector<const char*> charp_array(_it begin, _it end) {
	std::vector<const char*> v;
	v.reserve(end - begin);
	for (auto i = begin; i != end; ++i)
		v.push_back(i->c_str());
	return v;
}

#define MAX_LEN(A) max_len(A, sizeof(A)/sizeof(A[0]))

std::string convert_size(size_t size);

namespace Util { namespace String {

// Workaround since sprintf is inconsistent in double rounding for different implementations.
inline int format_double(double x, char *p, int64_t buf_size) {
	if (x >= 100.0)
		return snprintf(p, buf_size, "%lli", (long long)std::floor(x)); // for keeping output compatible with BLAST
	long long i = std::llround(x*10.0);
	return snprintf(p, buf_size, "%lli.%lli", i / 10, i % 10);
}

std::string replace(const std::string& s, char a, char b);
std::string ratio_percentage(const double x, const double y);
std::string ratio_percentage(const size_t x, const size_t y);
int64_t interpret_number(const std::string& s);
std::vector<std::string> tokenize(const char* str, const char* delimiters);
std::set<int32_t> parse_csv(const std::string& s);

template<typename T>
T convert_string(const char* s);
template<typename T>
T convert_string(const std::string& s) {
	return convert_string<T>(s.c_str());
}

std::string format(double number);

}}
