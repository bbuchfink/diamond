/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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
