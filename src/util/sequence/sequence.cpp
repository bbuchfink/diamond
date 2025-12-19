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

#include <stdexcept>
#include "sequence.h"
#include "stats/score_matrix.h"
#include "translate.h"
#include "../table.h"
#include "../string/tokenizer.h"

#define ARRAY_SIZE(x) sizeof(x)/sizeof(x[0])

using std::vector;
using std::array;
using std::string;
using std::min;
using std::runtime_error;

namespace Util { namespace Seq {

const char* id_delimiters = " \a\b\f\n\r\t\v\1";
const char* FASTA_HEADER_SEP[2] = {"\1", " >"};

void format(Sequence seq, const char *id, const char *qual, TextBuffer &out, const std::string &format, const ValueTraits& value_traits, Loc wrap) {
	const size_t l = seq.length();
	if(format == "fasta") {
		out << '>' << id << '\n';
		for (Loc i = 0; i < seq.length(); i += wrap) {
			seq.print(out, i, std::min(i + wrap, seq.length()), value_traits);
			out << '\n';
		}
	}
	else if (format == "fastq") {
		out << '@' << id << '\n';
		seq.print(out, value_traits);
		out << "\n+\n";
		out << qual << '\n';
	}
	else
		throw runtime_error("Invalid sequence file format");
}

std::string all_seqids(const char* s)
{
	string r, t;
	Util::String::Tokenizer<Util::String::StringDelimiters> tok(s, Util::String::StringDelimiters(FASTA_HEADER_SEP, ARRAY_SIZE(FASTA_HEADER_SEP)));
	while(tok.good()) {
		if (tok.ptr() != s)
			r.append("\1");
		tok >> t;
		r.append(seqid(t.c_str()));
	}
	return r;
}

string seqid(const char* title) {
	return string(title, title + find_first_of(title, id_delimiters));
}

string get_accession(const string& title, AccessionParsing& stat) {
	size_t i;
	string t(title);
	if (t.compare(0, 6, "UniRef") == 0) {
		t.erase(0, t.find('_', 0) + 1);
		++stat.uniref_prefix;
	}
	else if ((i = t.find_first_of('|', 0)) != string::npos) {
		if (t.compare(0, 3, "gi|") == 0) {
			t.erase(0, t.find_first_of('|', i + 1) + 1);
			i = t.find_first_of('|', 0);
			++stat.gi_prefix;
		}
		t.erase(0, i + 1);
		++stat.prefix_before_pipe;
		i = t.find_first_of('|', 0);
		if (i != string::npos) {
			t.erase(i);
			++stat.suffix_after_pipe;
		}
	}
	i = t.find_last_of('.');
	if (i != string::npos) {
		t.erase(i);
		++stat.suffix_after_dot;
	}
	return t;
}

vector<string> accession_from_title(const char* title, bool parse_seqids, AccessionParsing& stat)
{
	Util::String::Tokenizer<Util::String::StringDelimiters> tok(title, Util::String::StringDelimiters(FASTA_HEADER_SEP, ARRAY_SIZE(FASTA_HEADER_SEP)));
	string s;
	vector<string> r;
	while (tok.good()) {
		tok >> s;
		r.push_back(parse_seqids ? get_accession(Util::Seq::seqid(s.c_str()), stat) : Util::Seq::seqid(s.c_str()));
	}
	return r;
}

std::ostream& operator<<(std::ostream& s, const AccessionParsing& stat) {
	Util::Table t;
	t("UniRef prefix", stat.uniref_prefix);
	t("gi|xxx| prefix", stat.gi_prefix);
	t("xxx| prefix", stat.prefix_before_pipe);
	t("|xxx suffix", stat.suffix_after_pipe);
	t(".xxx suffix", stat.suffix_after_dot);
	t(":PDB= suffix", stat.pdb_suffix);
	s << t;
	return s;
}

static const char* TAB_ERR = "Tabulator character in sequence title";
static const char* SPACES_ERR = "Leading spaces in sequence title";
static const char* BLANK_ERR = "Blank sequence title";

const char* fix_title(string& s) {
	size_t i = 0;
	const char* r = nullptr;
	while (i < s.length() && s[i] < 33) ++i;
	if (i > 0) {
		s.erase(0, i);
		r = SPACES_ERR;
	}
	if (s.empty()) {
		s = "N/A";
		return BLANK_ERR;
	}
	for (size_t i = 0; i < s.length(); ++i) {
		if (s[i] == '\t') {
			s.replace(i, 1, "\\t");
			r = TAB_ERR;
		}
	}
	return r;
}

void get_title_def(const std::string& s, std::string& title, std::string& def)
{
	const size_t i = find_first_of(s.c_str(), id_delimiters);
	title = s.substr(0, i);
	if (i >= s.length())
		def.clear();
	else
		def = s.substr(i + 1);
}

bool is_fully_masked(const Sequence& seq) {
	Loc n = 0;
	for (const Letter* p = seq.data(); p < seq.end(); ++p)
		if (*p >= TRUE_AA)
			++n;
	return n == seq.length();
}

array<vector<Letter>, 6> translate(const Sequence& seq) {
	array<vector<Letter>, 6> out;
	if (seq.length() < 3)
		return out;
	Translator::translate(seq, out.data());
	return out;
}

Loc find_orfs(vector<Letter>& seq, const Loc min_len) {
	vector<Letter>::iterator it, begin = seq.begin();
	Loc n = 0;
	while ((it = std::find(begin, seq.end(), STOP_LETTER)) != seq.end()) {
		const Loc l = Loc(it - begin);
		if (l < min_len)
			std::fill(begin, it, MASK_LETTER);
		else
			n += l;
		begin = it + 1;
	}
	const Loc l = Loc(seq.end() - begin);
	if (l < min_len)
		std::fill(begin, seq.end(), MASK_LETTER);
	else
		n += l;
	return n;
}

bool looks_like_dna(const Sequence& seq) {
	array<Loc, AMINO_ACID_COUNT> count;
	count.fill(0);
	for (Loc i = 0; i < seq.length(); ++i)
		++count[(int)seq[i]];
	return count[(int)value_traits.from_char('A')]
		+ count[(int)value_traits.from_char('C')]
		+ count[(int)value_traits.from_char('G')]
		+ count[(int)value_traits.from_char('T')]
		+ count[(int)value_traits.from_char('N')] == seq.length();
}

std::vector<Score> window_scores(Sequence seq1, Sequence seq2, Loc window) {
	assert(seq1.length() == seq2.length());
	vector<Score> v;
	v.reserve(seq1.length());
	Score s = 0;
	const Loc l = min(seq1.length(), window);
	for (Loc i = 0; i < l; ++i) {
		s += score_matrix(seq1[i], seq2[i]);
		v.push_back(s);
	}
	Loc j = 0;
	for (Loc i = window; i < seq1.length(); ++i, ++j) {
		s += score_matrix(seq1[i], seq2[i]);
		s -= score_matrix(seq1[j], seq2[j]);
		v.push_back(s);
	}
	return v;
}

std::vector<Letter> from_string(const char* s, const ValueTraits& t, int64_t line) {
	vector<Letter> v;
	v.reserve(strlen(s));
	while (*s) {
		v.push_back(t.from_char(*s));
		++s;
	}
	return v;
}

void from_string(const char* s, vector<Letter>& out, const ValueTraits& t, int64_t line) {
	out.clear();
	out.reserve(strlen(s));
	while (*s) {
		out.push_back(t.from_char(*s));
		++s;
	}
}

void from_string(const string& s, vector<Letter>& out, const ValueTraits& t, int64_t line) {
	out.clear();
	out.reserve(s.length());
	for (char c : s) {
		if (c == '\n' || c == '\r')
			continue;
		out.push_back(t.from_char(c));
	}
}

std::string remove_newlines(const std::string& s) {
	std::string r;
	r.reserve(s.length());
	for (char c : s) {
		if (c != '\n' && c != '\r')
			r.push_back(c);
	}
	return r;
}

}}