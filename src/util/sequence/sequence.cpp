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

#include <algorithm>
#include <stdexcept>
#include "sequence.h"
#include "../util.h"
#include "../../stats/score_matrix.h"
#include "translate.h"

using std::vector;
using std::array;
using std::string;
using std::min;

namespace Util { namespace Seq {

const char* id_delimiters = " \a\b\f\n\r\t\v\1";

void format(Sequence seq, const char *id, const char *qual, OutputFile &out, const std::string &format, const ValueTraits& value_traits) {
	static TextBuffer buf;
	constexpr Loc WRAP = 160;
	const size_t l = seq.length();
	if(format == "fasta") {
		buf << '>' << id << '\n';
		for (Loc i = 0; i < seq.length(); i += WRAP) {
			seq.print(buf, i, std::min(i + WRAP, seq.length()), value_traits);
			buf << '\n';
		}
	}
	else if (format == "fastq") {
		buf << '@' << id << '\n';
		seq.print(buf, value_traits);
		buf << "\n+\n";
		buf << qual << '\n';
	}
	else
		throw std::runtime_error("Invalid sequence file format");
	out.write(buf.data(), buf.size());
	buf.clear();
}

std::string all_seqids(const char* s)
{
	string r;
	const vector<string> t(tokenize(s, "\1"));
	for (vector<string>::const_iterator i = t.begin(); i != t.end(); ++i) {
		if (i != t.begin())
			r.append("\1");
		r.append(i->substr(0, find_first_of(i->c_str(), id_delimiters)));
	}
	return r;
}

std::string seqid(const char* title, bool short_seqids)
{
	string s(title, title + find_first_of(title, id_delimiters));
	if (!short_seqids)
		return s;

	size_t i = s.find_first_of('|', 0);
	if (i != string::npos) {
		s.erase(0, i + 1);
		i = s.find_first_of('|', 0);
		if (i != string::npos) {
			if (i == 0)
				s.erase(0, 1);
			else
				s.erase(i);
		}
		return s;
	}
	else
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

}}