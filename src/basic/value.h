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
enum class SequenceType : int32_t { amino_acid=0, nucleotide=1 };
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

struct ValueTraits
{
	ValueTraits(const char *alphabet, Letter mask_char, const char *ignore, SequenceType seq_type);
	const char *alphabet;
	unsigned alphabet_size;
	Letter mask_char;
	Char_representation from_char;
	SequenceType seq_type;
};

#define AMINO_ACID_ALPHABET "ARNDCQEGHILKMFPSTWYVBJZX*_"
#define AMINO_ACID_COUNT (int(sizeof(AMINO_ACID_ALPHABET) - 1))

#define NUCLEOTIDE_ALPHABET "ACGTN"
#define NUCLEOTIDE_COUNT (sizeof(NUCLEOTIDE_ALPHABET) -1)

constexpr Letter MASK_LETTER = 23;
constexpr Letter STOP_LETTER = 24;
constexpr Letter SUPER_HARD_MASK = 25;
constexpr Letter DELIMITER_LETTER = 31;
constexpr Letter LETTER_MASK = 31;
constexpr Letter SEED_MASK = -128;
constexpr int32_t TRUE_AA = 20;

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
#ifdef SEQ_MASK
	return _mm_and_si128(x, _mm_set1_epi8(LETTER_MASK));
#else
	return x;
#endif
}
#endif

#ifdef __AVX2__
static inline __m256i letter_mask(__m256i x) {
#ifdef SEQ_MASK
	return _mm256_and_si256(x, _mm256_set1_epi8(LETTER_MASK));
#else
	return x;
#endif
}
#endif

extern const ValueTraits amino_acid_traits;
extern const ValueTraits nucleotide_traits;
extern ValueTraits value_traits;
extern ValueTraits input_value_traits;

static inline char to_char(Letter a)
{
	return value_traits.alphabet[(long)a];
}

struct AlignMode
{
	AlignMode(unsigned mode);
	static unsigned from_command(unsigned command);
	int check_context(int i) const
	{
		if (i >= query_contexts)
			throw std::runtime_error("Sequence context is out of bounds.");
		return i;
	}
	enum { blastp = 2, blastx = 3, blastn = 4 };
	SequenceType sequence_type, input_sequence_type;
	int mode, query_contexts, query_len_factor;
	bool query_translated;
};

extern AlignMode align_mode;
extern const Letter IUPACAA_TO_STD[32];
extern const Letter NCBI_TO_STD[28];

enum class Alphabet { STD, NCBI };

template<typename It>
void alph_ncbi_to_std(const It begin, const It end) {
	for (It i = begin; i != end; ++i) {
		const size_t l = *i;
		if (l >= sizeof(NCBI_TO_STD))
			throw std::runtime_error("Unrecognized sequence character in BLAST database");
		*i = NCBI_TO_STD[l];
	}
}

using Loc = int32_t;
using BlockId = int32_t;
using OId = int64_t;
using DictId = int64_t;
using Score = int32_t;
using TaxId = int32_t;
using CentroidId = OId;
using SuperBlockId = int32_t;