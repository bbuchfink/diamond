#ifndef UTIL_STRING_STRING_H_
#define UTIL_STRING_STRING_H_

#include <string>
#include <string.h>
#include <algorithm>
#include <ostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <stdio.h>

inline bool ends_with(const std::string &s, const char *t) {
	const size_t l = strlen(t);
	return s.length() >= l && strncmp(&s[s.length() - l], t, l) == 0;
}

inline std::string& rstrip(std::string &s, const char *t) {
	const size_t l = strlen(t);
	if (s.length() < l)
		return s;
	if (s.compare(s.length() - l, std::string::npos, t) == 0)
		return s.erase(s.length() - l, std::string::npos);
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
inline int format_double(double x, char *p) {
	if (x >= 100.0)
		return sprintf(p, "%lli", (long long)std::floor(x)); // for keeping output compatible with BLAST
	long long i = std::llround(x*10.0);
	return sprintf(p, "%lli.%lli", i / 10, i % 10);
}


}}

#endif