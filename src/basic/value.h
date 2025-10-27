/****
Copyright � 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

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
#include <stdexcept>
#include <string>
#include <cstdint>
#include "util/simd.h"

typedef signed char Letter;
enum class SequenceType : int32_t { amino_acid = 0, nucleotide = 1 };

struct CharRepresentation
{
	CharRepresentation(unsigned size, const char* chars, char mask, const char* mask_chars);
	Letter operator()(char c) const
	{
		if (data_[(long)c] == invalid)
			throw std::runtime_error("Invalid character in sequence: " + (c >= 32 && c < 127 ? std::string("'") + c + "'" : "ASCII " + std::to_string((long)c)));
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
	CharRepresentation from_char;
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

#ifdef __ARM_NEON
static inline int8x16_t letter_mask(int8x16_t x) {
#ifdef SEQ_MASK
	return vandq_s8(x, vdupq_n_s8(LETTER_MASK));
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
	const char* to_string() const {
		static const char* mode_str[] = { 0, 0, "blastp", "blastx", "blastn" };
		return mode_str[mode];
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
using OId = int_fast64_t;
using DictId = int64_t;
using Score = int32_t;
using TaxId = int32_t;
using CentroidId = OId;
using SuperBlockId = int32_t;