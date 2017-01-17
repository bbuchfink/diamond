/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#ifndef VALUE_H_
#define VALUE_H_

#include <assert.h>
#include <string.h>
#include <sstream>
#include <stdexcept>
#include "const.h"
#include "config.h"
#include "../util/util.h"

typedef char Letter;
typedef enum { amino_acid=0, nucleotide=1 } Sequence_type;
struct Amino_acid {};
struct Nucleotide {};

struct invalid_sequence_char_exception : public std::exception
{
	const std::string msg;
	invalid_sequence_char_exception(char ch) :
		msg(std::string("Invalid character (") + print_char(ch) + ") in sequence")
	{ }
	~invalid_sequence_char_exception() throw()
	{ }
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
};

struct Char_representation
{
	Char_representation(unsigned size, const char *chars, char mask, const char *mask_chars)
	{
		memset(data_, invalid, sizeof(data_));
		for (unsigned i = 0; i<size; ++i) {
			assert(chars[i] != (char)invalid);
			data_[(long)chars[i]] = i;
			data_[(long)tolower(chars[i])] = i;
		}
		while (*mask_chars != 0) {
			const char ch = *mask_chars;
			data_[(long)ch] = mask;
			data_[(long)tolower(ch)] = mask;
			++mask_chars;
		}
	}
	Letter operator()(char c) const
	{
		if (data_[(long)c] == invalid)
			throw invalid_sequence_char_exception(c);
		return data_[(long)c];
	}
private:
	static const char invalid;
	Letter data_[256];
};

struct Value_traits
{
	Value_traits(const char *alphabet, Letter mask_char, const char *ignore);	
	const char *alphabet;
	unsigned alphabet_size;
	Letter mask_char;
	Char_representation from_char;
};

extern const Value_traits amino_acid_traits;
extern const Value_traits nucleotide_traits;
extern Value_traits value_traits;
extern Value_traits input_value_traits;

inline char to_char(Letter a)
{
	return value_traits.alphabet[(long)a];
}

struct Align_mode
{
	Align_mode(unsigned mode);
	static unsigned from_command(unsigned command);
	unsigned check_context(unsigned i) const
	{
		if (i >= query_contexts)
			throw std::runtime_error("Sequence context is out of bounds.");
		return i;
	}
	enum { blastp = 2, blastx = 3, blastn = 4 };
	Sequence_type sequence_type, input_sequence_type;
	unsigned mode, query_contexts;
	int query_len_factor;
	bool query_translated;
};

extern Align_mode align_mode;
extern const double background_freq[20];

#endif /* VALUE_H_ */
