#ifndef UTIL_STRING_STRING_H_
#define UTIL_STRING_STRING_H_

#include <string>
#include <string.h>
#include <algorithm>
#include <ostream>

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

#define MAX_LEN(A) max_len(A, sizeof(A)/sizeof(A[0]))

std::string convert_size(size_t size);

#endif