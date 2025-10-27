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
#include "basic/sequence.h"
#include "string_set.h"
#include "basic/value.h"

struct SequenceSet : public StringSetBase<Letter, Sequence::DELIMITER, 1>
{

	SequenceSet(Alphabet alphabet = Alphabet::STD);
	SequenceSet(StringSetBase&& string_set):
		StringSetBase(string_set),
		alphabet_(Alphabet::STD)
	{}
	
	void print_stats() const;

	Sequence operator[](size_t i) const
	{
		return Sequence(ptr(i), (Loc)length(i));
	}

	std::pair<Length, Length> len_bounds(Length min_len) const;

	Length max_len(size_t begin, size_t end) const;

	std::vector<uint32_t> partition(unsigned n_part, bool shortened = false, bool context_reduced = false) const;

	size_t reverse_translated_len(size_t i) const;

	TranslatedSequence translated_seq(const Sequence& source, size_t i) const;

	size_t avg_len() const;

	virtual ~SequenceSet();

	Alphabet alphabet() const {
		return alphabet_;
	}

	Alphabet& alphabet() {
		return alphabet_;
	}

	void convert_to_std_alph(size_t id);
	void convert_all_to_std_alph(size_t threads);
	std::vector<std::pair<Loc, BlockId>> lengths() const;
	Loc source_length(BlockId i) const;
	
private:

	Alphabet alphabet_;

};

size_t max_id_len(const StringSet& ids);