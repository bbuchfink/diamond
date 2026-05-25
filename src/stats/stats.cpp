/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <algorithm>
#include <iterator>
#include "stats/score_matrix.h"
#include "standard_matrix.h"
#include "matrices/blosum45.h"
#include "matrices/blosum50.h"
#include "matrices/blosum62.h"
#include "matrices/blosum80.h"
#include "matrices/blosum90.h"
#include "matrices/pam30.h"
#include "matrices/pam70.h"
#include "matrices/pam250.h"
#include "basic/reduction.h"

using std::string;
using std::max;
using std::min;
using std::map;
using std::runtime_error;

const ValueTraits amino_acid_traits(AMINO_ACID_ALPHABET, 23, "UO-", SequenceType::amino_acid);
const ValueTraits nucleotide_traits("ACGTN", 4, "MRWSYKVHDBX", SequenceType::nucleotide);
ValueTraits value_traits(amino_acid_traits);
ValueTraits input_value_traits(amino_acid_traits);
Reduction Reduction::instance("A KR EDNQ C G H ILVM FYW P ST");

namespace Search {
	Reduction murphy10("A KR EDNQ C G H ILVM FYW P ST");
	Reduction steinegger12("AST C DN EQ FY G H IV KR LM P W");
	Reduction no_reduction("A S T C D N E Q F Y G H I V K R L M P W");
	Reduction dna("A C G T");
}

namespace Stats {

const map<string, const StandardMatrix&> StandardMatrix::matrices = { { "blosum45", blosum45 }, { "blosum62", blosum62 },
	{ "blosum50", blosum50 }, { "blosum80", blosum80 }, { "blosum90", blosum90 },
	{ "pam250", pam250 }, { "pam30", pam30 }, { "pam70", pam70 } };

const StandardMatrix& StandardMatrix::get(const string& name) {
	string n;
	std::transform(name.begin(), name.end(), std::back_inserter(n), [](unsigned char c) { return tolower(c); });
	auto it = matrices.find(n);
	if (it == matrices.end())
		throw std::runtime_error("Unknown scoring matrix: " + name);
	return it->second;
}

const StandardMatrix::Parameters& StandardMatrix::constants(int gap_exist, int gap_extend) const {
	const double g = gap_exist, e = gap_extend;
	for (const auto& i : parameters)
		if (i.gap_exist == g && i.gap_extend == e)
			return i;
	throw runtime_error("Gap penalty settings are outside the supported range for this scoring matrix.");
}

const StandardMatrix::Parameters& StandardMatrix::ungapped_constants() const {
	return parameters.front();
}

/*double approx_id(int score, Loc range1, Loc range2) {
	static const std::array<std::pair<double, double>, 20> lookup = { {
		{ 80.0, 1.555513051244674},
		{ 81.0, 1.572942897590095},
		{ 82.0, 1.58750168603496},
		{ 83.0, 1.6000124264834454},
		{ 84.0, 1.6161145696530383 },
		{ 85.0, 1.6372885235111332 },
		{ 86.0, 1.6580185244345147 },
		{ 87.0, 1.6790883718225538 },
		{ 88.0, 1.700225363839303 },
		{ 89.0, 1.7246740739439146 },
		{ 90.0, 1.7479513469270396 },
		{ 91.0, 1.7727672853690373 },
		{92.0, 1.792468675471071},
		{93.0, 1.8155660510729361},
		{94.0, 1.8343271923898516},
		{95.0, 1.8594309345637143},
		{96.0,1.8803024476343855},
	{ 97.0,1.9063177406803247 },
	{ 98.0, 1.9286849708679306 },
		{99.0, 1.9560011599829659} } };
	const Loc m = max(range1, range2);
	if (m == 0)
		return 100.0;
	const double b = bit_score / m;
	for (auto i = lookup.rbegin(); i != lookup.rend(); ++i)
		if (b >= i->second)
			return i->first;
	return 0.0;
}*/

double approx_id(Score raw_score, Loc range1, Loc range2) {
	const Loc m = max(range1, range2);
	if (m == 0)
		return 100.0;
	return min(max((double)raw_score / m * 16.56 + 11.41, 0.0), 100.0);
}

double approx_id(int score, Interval query_range, Interval target_range, const Sequence& query, const Sequence& target) {
	
	if (query.length() == 0 || target.length() == 0)
		throw runtime_error("Cannot compute approximate identity for empty sequences.");
	Score s = 0;
	if (query_range.length() < target_range.length()) {
		for (int i = query_range.begin_; i < query_range.end_; ++i) {
			const Letter l = query[i];
			s += score_matrix(l, l);
		}
		return (double)score / s * 100.0;
	}
	else {
		for (int i = target_range.begin_; i < target_range.end_; ++i) {
			const Letter l = target[i];
			s += score_matrix(l, l);
		}
		return (double)score / s * 100.0;
	}
}

}