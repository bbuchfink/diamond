#ifndef UTIL_STRING_TOKENIZER_H_
#define UTIL_STRING_TOKENIZER_H_

#include <string>
#include <string.h>
#include <stdexcept>
#include <stdlib.h>

namespace Util { namespace String {

struct TokenizerException : public std::runtime_error {
	TokenizerException():
		std::runtime_error("Tokenizer Exception")
	{}
	TokenizerException(const std::string &msg):
		std::runtime_error(msg)
	{}
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
		if (!good())
			throw TokenizerException();
		char *end;
		x = strtol(p, &end, 10);
		if (end == p)
			throw TokenizerException();
		if (strncmp(end, delimiter, len) == 0)
			p = end + len;
		else {
			if (*end != '\0')
				throw TokenizerException();
			p = nullptr;
		}
		return *this;
	}

	Tokenizer& operator>>(int &x) {
		long y;
		*this >> y;
		x = (int)y;
		return *this;
	}

	Tokenizer& operator>>(double &x) {
		if (!good())
			throw TokenizerException("No token left");
		char *end;
		x = strtod(p, &end);
		if (end == p)
			throw TokenizerException("Unable to parse double");
		if (strncmp(end, delimiter, len) == 0)
			p = end + len;
		else {
			if (*end != '\0')
				throw TokenizerException("Invalid char in double");
			p = nullptr;
		}
		return *this;
	}

	Tokenizer& operator>>(float &x) {
		if (!good())
			throw TokenizerException("No token left");
		char *end;
		x = strtof(p, &end);
		if (end == p)
			throw TokenizerException("Unable to parse float");
		if (strncmp(end, delimiter, len) == 0)
			p = end + len;
		else {
			if (*end != '\0')
				throw TokenizerException("Invalid char in float");
			p = nullptr;
		}
		return *this;
	}

	bool good() const {
		return p != nullptr && *p != '\0';
	}

	void skip_to(char c) {
		const char* q = strchr(p, c);
		p = q ? q + 1 : nullptr;
	}

private:
	const char *p, *delimiter;
	size_t len;

};

}}

#endif
