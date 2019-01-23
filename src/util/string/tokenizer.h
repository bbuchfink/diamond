#ifndef UTIL_STRING_TOKENIZER_H_
#define UTIL_STRING_TOKENIZER_H_

#include <string>
#include <string.h>
#include <stdexcept>
#include <stdlib.h>

namespace Util { namespace String {

struct TokenizerException : public std::exception {
};

struct Skip {};

struct Tokenizer {

	Tokenizer(const std::string &s, const char *delimiter):
		p(s.c_str()),
		delimiter(delimiter),
		len(strlen(delimiter))
	{}

	Tokenizer& operator>>(const Skip&) {
		if (p == nullptr)
			throw TokenizerException();
		const char *d = strstr(p, delimiter);
		if (d) {
			p = d + len;
		}
		else {
			p = nullptr;
		}
		return *this;
	}

	Tokenizer& operator>>(std::string &s) {
		if (p == nullptr)
			throw TokenizerException();
		const char *d = strstr(p, delimiter);
		if (d) {
			s.assign(p, d - p);
			p = d + len;
		}
		else {
			s.assign(p);
			p = nullptr;
		}
		return *this;
	}

	Tokenizer& operator>>(long int &x) {
		if (p == nullptr || *p == '\0')
			throw TokenizerException();
		char *end;
		x = strtol(p, &end, 10);
		if (end == p)
			throw TokenizerException();
		if (strncmp(end, delimiter, len) == 0)
			p = end + len;
		else {
			if (*p != '\0')
				throw TokenizerException();
			p = nullptr;
		}
		return *this;
	}

private:
	const char *p, *delimiter;
	size_t len;

};

}}

#endif
