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

#include <cassert>
#include <string.h>
#include "value.h"

const Letter CharRepresentation::invalid = '\xff';

CharRepresentation::CharRepresentation(unsigned size, const char* chars, char mask, const char* mask_chars)
{
	memset(data_, invalid, sizeof(data_));
	for (unsigned i = 0; i < size; ++i) {
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

ValueTraits::ValueTraits(const char* alphabet, Letter mask_char, const char* ignore, SequenceType seq_type) :
	alphabet(alphabet),
	alphabet_size((unsigned)strlen(alphabet)),
	mask_char(mask_char),
	from_char(CharRepresentation((unsigned)alphabet_size, alphabet, mask_char, ignore)),
	seq_type(seq_type)
{}

// 15 = O, 21 = U
const Letter IUPACAA_TO_STD[32] = { -1, 0, 20, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10, 12, 2, MASK_LETTER, 14, 5, 1, 15, 16, MASK_LETTER, 19, 17, 23, 18, 22, -1, -1, -1, -1, 24 };
// 24 = U, 26 = O
const Letter NCBI_TO_STD[28] = { MASK_LETTER, 0, 20, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 23, 18, 22, MASK_LETTER, 24, MASK_LETTER, 21 };