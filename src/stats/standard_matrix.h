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

#pragma once

#include <map>
#include <array>
#include <limits.h>
#include <vector>
#include <float.h>
#include "../basic/value.h"

namespace Stats {

const double INT2_MAX = DBL_MAX;
const size_t NCBI_ALPH = 28;
using FreqRatios = double[NCBI_ALPH][NCBI_ALPH];

struct StandardMatrix {

	struct Parameters {
		double gap_exist;
		double gap_extend;
		double reserved;
		double Lambda;
		double K;
		double H;
		double alpha;
		double beta;
		double C;
		double alpha_v;
		double sigma;
	};

	int default_gap_exist, default_gap_extend;
	std::vector<Parameters> parameters;
	std::array<int8_t, AMINO_ACID_COUNT * AMINO_ACID_COUNT> scores;
	double joint_probs[TRUE_AA][TRUE_AA];
	std::array<double, TRUE_AA> background_freqs;
	FreqRatios freq_ratios;

	const Parameters& constants(int gap_exist, int gap_extend) const;
	const Parameters& ungapped_constants() const;
	static const StandardMatrix& get(const std::string& name);
	static const std::map<std::string, const StandardMatrix&> matrices;

};

extern const StandardMatrix blosum45, blosum50, blosum62, blosum80, blosum90, pam250, pam30, pam70;

}