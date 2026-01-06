/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
						Benjamin Buchfink

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include <string>
#include <string.h>
#include <stdexcept>
#include <stdlib.h>
#include <tuple>

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

struct CharDelimiter {
	CharDelimiter (char c):
		c(c)
	{}
	std::pair<const char*, const char*> scan(const char* p) const {
		const char* q = strchr(p, c);
		return q ? std::pair<const char*, const char*>(q, q + 1) : std::pair<const char*, const char*>(nullptr, nullptr);
	}
	const char* next(const char* p) const {
		if (p[0] == c)
			return p + 1;
		else {
			if (*p != '\0')
				throw TokenizerException();
			return nullptr;
		}
	}
	const char c;
};

struct StringDelimiter {
	StringDelimiter(const char* s):
		s(s),
		len((int)strlen(s))
	{}
	std::pair<const char*, const char*> scan(const char* p) const {
		const char* q = strstr(p, s);
		return q ? std::pair<const char*, const char*>(q, q + len) : std::pair<const char*, const char*>(nullptr, nullptr);
	}
	const char* next(const char* p) const {
		if (strncmp(p, s, len) == 0)
			return p + len;
		else {
			if (*p != '\0')
				throw TokenizerException();
			return nullptr;
		}
	}
	const char* s;
	const int len;
};

struct StringDelimiters {
	StringDelimiters(const char** s, int n) :
		s(s),
		n(n)
	{}
	std::pair<const char*, const char*> scan(const char* p) const {
		for (int i = 0; i < n; ++i) {
			const char* ptr = strstr(p, s[i]);
			if (ptr)
				return { ptr,ptr + strlen(s[i]) };
		}
		return { nullptr,nullptr };
	}
	const char* next(const char* p) const {
		for (int i = 0; i < n; ++i) {
			const size_t l = strlen(s[i]);
			if (strncmp(p, s[i], l) == 0)
				return p + l;
		}
		if (*p != '\0')
			throw TokenizerException();
		return nullptr;
	}
	const char** s;
	const int n;
};

template<typename Delimiters>
struct Tokenizer {

	Tokenizer(const std::string &s, Delimiters delimiters):
		p(s.c_str()),
		delimiters(delimiters)
	{}

	Tokenizer(const char* s, Delimiters delimiters):
		p(s),
		delimiters(delimiters)
	{}

	const char* ptr() const {
		return p;
	}

	void set(const char* ptr) {
		p = ptr;
	}

	Tokenizer& operator>>(const Skip&) {
		if (p == nullptr)
			throw TokenizerException();
		std::tie(std::ignore, p) = delimiters.scan(p);
		return *this;
	}

	Tokenizer& operator>>(std::string &s) {
		if (p == nullptr)
			throw TokenizerException();
		const auto d = delimiters.scan(p);
		if (d.first)
			s.assign(p, d.first - p);
		else
			s.assign(p);
		p = d.second;
		return *this;
	}

	Tokenizer& operator>>(int64_t &x) {
		if (!good())
			throw TokenizerException();
		char *end;
		const long long n = strtoll(p, &end, 10);
		x = (int64_t)n;
		if (end == p)
			throw TokenizerException();
		p = delimiters.next(end);
		return *this;
	}
	
	Tokenizer& operator>>(uint64_t& x) {
		if (!good())
			throw TokenizerException();
		char* end;
		const long long n = strtoll(p, &end, 10);
		x = (uint64_t)n;
		if (end == p)
			throw TokenizerException();
		p = delimiters.next(end);
		return *this;
	}

	Tokenizer& operator>>(int &x) {
		int64_t y;
		*this >> y;
		x = (int)y;
		return *this;
	}

	Tokenizer& operator>>(unsigned& x) {
		int64_t y;
		*this >> y;
		x = (unsigned)y;
		return *this;
	}

	Tokenizer& operator>>(double &x) {
		if (!good())
			throw TokenizerException("No token left");
		char *end;
		x = strtod(p, &end);
		if (end == p)
			throw TokenizerException("Unable to parse double");
		p = delimiters.next(end);
		return *this;
	}

	Tokenizer& operator>>(float &x) {
		if (!good())
			throw TokenizerException("No token left");
		char *end;
		x = strtof(p, &end);
		if (end == p)
			throw TokenizerException("Unable to parse float");
		p = delimiters.next(end);
		return *this;
	}

	bool good() const {
		return p != nullptr && *p != '\0';
	}

	void skip_to(char c) {
		const char* q = strchr(p, c);
		p = q ? q + 1 : nullptr;
	}

	std::string getline() const {
		const char* q = strchr(p, '\n');
		if (q == nullptr)
			q = strchr(p, '\0');
		return std::string(p, q);
	}

private:

	const char* p;
	const Delimiters delimiters;

};

}}