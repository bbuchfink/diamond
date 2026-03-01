/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <string>
#include <array>
#include <vector>
#include "basic/sequence.h"

namespace Util { namespace Seq {

void format(Sequence seq, const char *id, const char *qual, TextBuffer &out, const std::string &format, const ValueTraits& value_traits, Loc wrap = 160);

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
std::string remove_newlines(const std::string& s);

}}

template<typename Callback>
void read_fasta(std::istream& in, Callback& f)
{
	std::string line, id, seq;
	std::streampos pos = in.tellg(), start = pos;
	std::getline(in, line);
	if (line.empty() || line.front() != '>')
		throw std::runtime_error("FASTA format error: file does not start with '>'");
	do {
		if (!line.empty() && line.back() == '\r')
			line.pop_back();
		if (line.empty())
			continue;
		if (line[0] == '>') {			
			if (!id.empty())
				f(id, seq, start);
			start = pos;
			id = line.substr(1);
			if (id.empty())
				throw std::runtime_error("FASTA format error: empty id at file offset " + std::to_string(pos));
			seq.clear();
		}
		else {
			seq += line;
		}
		pos = in.tellg();
	} while (std::getline(in, line));
	if (!id.empty())
		f(id, seq, start);
}