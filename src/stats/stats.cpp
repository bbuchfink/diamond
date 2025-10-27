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

#include <algorithm>
#include <iterator>
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

double approx_id(Score raw_score, Loc range1, Loc range2) {
	const Loc m = max(range1, range2);
	if (m == 0)
		return 100.0;
	return min(max((double)raw_score / m * 16.56 + 11.41, 0.0), 100.0);
}

}