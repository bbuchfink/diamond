/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef VALUE_H_
#define VALUE_H_

#include "value_type.h"
#include "const.h"

typedef enum { amino_acid, nucleotide } Sequence_type;

Sequence_type input_sequence_type()
{ return program_options::command == program_options::blastp ? amino_acid : nucleotide; }

Sequence_type sequence_type()
{ return program_options::command == program_options::blastn ? nucleotide : amino_acid; }

size_t query_contexts()
{
	switch(program_options::command) {
	case program_options::blastn: return 2;
	case program_options::blastx: return 6;
	default: return 1;
	}
}

bool query_translated()
{ return program_options::command == program_options::blastx ? true : false; }

int query_len_factor()
{ return program_options::command == program_options::blastx ? 3 : 1; }

template<typename _val>
struct Char_representation
{
	Char_representation(unsigned size, const char *chars, char mask, const char *mask_chars)
	{
		memset(data_, invalid, sizeof(data_));
		for(unsigned i=0;i<size;++i) {
			assert(chars[i] != (char)invalid);
			data_[(long)chars[i]] = i;
			data_[(long)tolower(chars[i])] = i;
		}
		while(*mask_chars != 0) {
			const char ch = *mask_chars;
			data_[(long)ch] = mask;
			data_[(long)tolower(ch)] = mask;
			++mask_chars;
		}
	}
	_val operator()(char c) const
	{
		if(data_[(long)c] == invalid)
			throw invalid_sequence_char_exception (c);
		return data_[(long)c];
	}
private:
	static const _val invalid;
	_val data_[256];
};

template<> const Amino_acid Char_representation<Amino_acid>::invalid = 0xff;
template<> const Nucleotide Char_representation<Nucleotide>::invalid = 0xff;

template<typename _val>
struct Value_traits
{ };

template<>
struct Value_traits<Amino_acid>
{
	enum { ALPHABET_SIZE = 25 };
	static const Amino_acid				MASK_CHAR;
	static const char*					ALPHABET;
	static const Char_representation<Amino_acid>	from_char;
};

const Amino_acid					Value_traits<Amino_acid>::MASK_CHAR = 23;
const char* Value_traits<Amino_acid>::ALPHABET = "ARNDCQEGHILKMFPSTWYVBJZX*";
const Char_representation<Amino_acid> Value_traits<Amino_acid>::from_char (Value_traits<Amino_acid>::ALPHABET_SIZE, Value_traits<Amino_acid>::ALPHABET, Value_traits<Amino_acid>::MASK_CHAR, "UO");

template<>
struct Value_traits<const Amino_acid> : public Value_traits<Amino_acid>
{ };

template<>
struct Value_traits<Nucleotide>
{
	enum { ALPHABET_SIZE = 5 };
	static const Nucleotide				MASK_CHAR;
	static const char*					ALPHABET;
	static const Char_representation<Nucleotide>	from_char;
};

const Nucleotide Value_traits<Nucleotide>::MASK_CHAR = 4;
const char* Value_traits<Nucleotide>::ALPHABET = "ACGTN";
const Char_representation<Nucleotide> Value_traits<Nucleotide>::from_char (Value_traits<Nucleotide>::ALPHABET_SIZE, Value_traits<Nucleotide>::ALPHABET, Value_traits<Nucleotide>::MASK_CHAR, "MRWSYKVHDBX");

template<>
struct Value_traits<const Nucleotide> : public Value_traits<Nucleotide>
{ };

char to_char(Amino_acid a)
{ return Value_traits<Amino_acid>::ALPHABET[a]; }

template<>
struct Value_traits<char>
{
	static char from_char(char c)
	{ return c; }
};

#endif /* VALUE_H_ */
