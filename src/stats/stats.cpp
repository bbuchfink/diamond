/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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
#include <iterator>
#include <ctype.h>
#include "standard_matrix.h"

using std::string;
using std::max;
using std::min;

namespace Stats {

const std::map<std::string, const StandardMatrix&> StandardMatrix::matrices = { { "blosum45", blosum45 }, { "blosum62", blosum62 }, { "blosum50", blosum50 }, { "blosum80", blosum80 }, { "blosum90", blosum90 },
	{ "pam250", pam250 }, { "pam30", pam30 }, { "pam70", pam70 } };

const StandardMatrix& StandardMatrix::get(const std::string& name) {
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
	throw std::runtime_error("Gap penalty settings are not supported for this scoring matrix.");
}

const StandardMatrix::Parameters& StandardMatrix::ungapped_constants() const {
	return parameters.front();
}

double approx_id(Score raw_score, Loc range1, Loc range2) {
	return min(max((double)raw_score / max(range1, range2) * 16.56 + 11.41, 0.0), 100.0);
}

}