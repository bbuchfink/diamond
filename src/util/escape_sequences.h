/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef ESCAPE_SEQUENCES_H_
#define ESCAPE_SEQUENCES_H_

#include <limits.h>

struct EscapeSequence
{
	char c;
	const char *seq;
};

struct EscapeSequences
{

	EscapeSequences(const EscapeSequence seqs[], size_t n)
	{
		for (unsigned char c = 0; c < UCHAR_MAX; ++c) {
			seq_[(size_t)c][0] = c;
			seq_[(size_t)c][1] = '\0';
		}
		for (size_t i = 0; i < n; ++i)
			strcpy(seq_[(size_t)seqs[i].c], seqs[i].seq);
	}

	const char* escape(char c) const
	{
		return (char*)seq_[(size_t)c];
	}

	void escape(const char *s, size_t len, string &out) const
	{
		out.reserve(len);
		for (size_t i = 0; i < len; ++i)
			out += escape(*(s++));
	}

	void escape(const char *s, string &out) const
	{
		escape(s, strlen(s), out);
	}

	void escape(const string &s, string &out) const
	{
		escape(s.c_str(), s.length(), out);
	}

	static const EscapeSequences XML;

private:

	char seq_[UCHAR_MAX][7];

	static const EscapeSequence xml_data[5];

};

#endif