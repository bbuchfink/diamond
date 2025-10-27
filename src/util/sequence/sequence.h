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
#include <string>
#include <array>
#include <vector>
#include "basic/sequence.h"
#include "../io/output_file.h"

namespace Util { namespace Seq {

void format(Sequence seq, const char *id, const char *qual, TextBuffer &out, const std::string &format, const ValueTraits& value_traits);

static inline Sequence clip(const Letter *seq, int len, int anchor) {
	const Letter *a = seq + anchor, *begin = seq, *end = seq + len, *p;
	for(;;) {
		p = (const Letter*)memchr(begin, (int)Sequence::DELIMITER, end - begin);
		if (p == nullptr)
			return Sequence(begin, end);
		if (p >= a)
			return Sequence(begin, p);
		begin = p + 1;
	}
}

extern const char* id_delimiters;
extern const char* FASTA_HEADER_SEP[2];

struct AccessionParsing {
	AccessionParsing() :
		uniref_prefix(0),
		gi_prefix(0),
		prefix_before_pipe(0),
		suffix_after_pipe(0),
		suffix_after_dot(0),
		pdb_suffix(0)
	{}
	friend std::ostream& operator<<(std::ostream& s, const AccessionParsing& stat);
	int64_t uniref_prefix, gi_prefix, prefix_before_pipe, suffix_after_pipe, suffix_after_dot, pdb_suffix;
};

std::string get_accession(const std::string& t, AccessionParsing& stat);
std::vector<std::string> accession_from_title(const char* title, bool parse_seqids, AccessionParsing& stat);
std::string all_seqids(const char* s);
std::string seqid(const char* title);
std::vector<std::string> seq_titles(const char* title);

void get_title_def(const std::string& s, std::string& title, std::string& def);
bool is_fully_masked(const Sequence& seq);
std::array<std::vector<Letter>, 6> translate(const Sequence& seq);
Loc find_orfs(std::vector<Letter>& seq, const Loc min_len);
bool looks_like_dna(const Sequence& seq);
std::vector<Score> window_scores(Sequence seq1, Sequence seq2, Loc window);
const char* fix_title(std::string& s);
std::vector<Letter> from_string(const char* s, const ValueTraits& t, int64_t line);
void from_string(const char* s, std::vector<Letter>& out, const ValueTraits& t, int64_t line);
void from_string(const std::string& s, std::vector<Letter>& out, const ValueTraits& t, int64_t line);

}}