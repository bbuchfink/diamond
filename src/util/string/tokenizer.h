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

	const char* ptr() const {
		return p;
	}

	void set(const char* ptr) {
		p = ptr;
	}

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

	Tokenizer& operator>>(int64_t &x) {
		if (!good())
			throw TokenizerException();
		char *end;
		const long long n = strtoll(p, &end, 10);
		x = (int64_t)n;
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

	std::string getline() const {
		const char* q = strchr(p, '\n');
		if (q == nullptr)
			q = strchr(p, '\0');
		return std::string(p, q);
	}

private:
	const char *p, *delimiter;
	size_t len;

};

}}
