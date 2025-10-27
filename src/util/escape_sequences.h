/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "text_buffer.h"

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
		return (const char*)seq_[(size_t)c];
	}

	void escape(const char *s, size_t len, std::string &out) const
	{
		out.reserve(len);
		for (size_t i = 0; i < len; ++i)
			out += escape(*(s++));
	}

	void escape(const char *s, std::string &out) const
	{
		escape(s, strlen(s), out);
	}

	void escape(const std::string &s, std::string &out) const
	{
		escape(s.c_str(), s.length(), out);
	}

	static const EscapeSequences XML;

private:

	char seq_[UCHAR_MAX][7];

	static const EscapeSequence xml_data[5];

};

inline void print_escaped_until(TextBuffer& buf, const char* s, const char* delimiters, const EscapeSequences* esc)
{
    if (esc == 0)
        buf.write_until(s, delimiters);
    else {
        std::string tmp;
        esc->escape(s, find_first_of(s, delimiters), tmp);
        buf << tmp;
    }
}

inline void print_escaped(TextBuffer& buf, const std::string& s, const EscapeSequences* esc)
{
	if (esc == 0)
		buf << s;
	else {
		std::string tmp;
		esc->escape(s, tmp);
		buf << tmp;
	}
}