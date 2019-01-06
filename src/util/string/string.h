#ifndef UTIL_STRING_STRING_H_
#define UTIL_STRING_STRING_H_

#include <string>
#include <string.h>

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

#endif