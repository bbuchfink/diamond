/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include <stdexcept>
#include <string>
#include "const.h"
#include "../util/simd.h"

typedef signed char Letter;
typedef enum { amino_acid=0, nucleotide=1 } Sequence_type;
struct Amino_acid {};
struct Nucleotide {};

struct invalid_sequence_char_exception : public std::exception
{
	const std::string msg;
	invalid_sequence_char_exception(char ch);
	~invalid_sequence_char_exception() noexcept
	{ }
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
};

struct Char_representation
{
	Char_representation(unsigned size, const char* chars, char mask, const char* mask_chars);
	Letter operator()(char c) const
	{
		if (data_[(long)c] == invalid)
			throw invalid_sequence_char_exception(c);
		return data_[(long)c];
	}
private:
	static const Letter invalid;
	Letter data_[256];
};

struct Value_traits
{
	Value_traits(const char *alphabet, Letter mask_char, const char *ignore, Sequence_type seq_type);
	const char *alphabet;
	unsigned alphabet_size;
	Letter mask_char;
	Char_representation from_char;
	Sequence_type seq_type;
};

#define AMINO_ACID_ALPHABET "ARNDCQEGHILKMFPSTWYVBJZX*_"
#define AMINO_ACID_COUNT (sizeof(AMINO_ACID_ALPHABET) - 1)

constexpr Letter MASK_LETTER = 23;
constexpr Letter STOP_LETTER = 24;
constexpr Letter SUPER_HARD_MASK = 25;
constexpr Letter DELIMITER_LETTER = 31;
constexpr Letter LETTER_MASK = 31;
constexpr Letter SEED_MASK = -128;

static inline bool is_amino_acid(Letter x) {
	return x != MASK_LETTER && x != DELIMITER_LETTER && x != STOP_LETTER;
}

static inline Letter letter_mask(Letter x) {
#ifdef SEQ_MASK
	return Letter(x & LETTER_MASK);
#else
	return x;
#endif
}

#ifdef __SSE2__
static inline __m128i letter_mask(__m128i x) {
	return _mm_and_si128(x, SIMD::_mm_set1_epi8(LETTER_MASK));
}
#endif

#ifdef __AVX2__
static inline __m256i letter_mask(__m256i x) {
	return _mm256_and_si256(x, ::SIMD::_mm256_set1_epi8(LETTER_MASK));
}
#endif

extern const Value_traits amino_acid_traits;
extern const Value_traits nucleotide_traits;
extern Value_traits value_traits;
extern Value_traits input_value_traits;

static inline char to_char(Letter a)
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