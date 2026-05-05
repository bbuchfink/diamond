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

#include <vector>
#include <sstream>
#include "basic/config.h"
#include "util/string/string.h"

using std::stringstream;
using std::string;
using std::vector;
using std::runtime_error;

namespace Cluster {

vector<string> cluster_steps(double approx_id, bool linear) {
	if (!config.cluster_steps.empty()) {
		for (vector<string>::const_iterator it = config.cluster_steps.begin() + 1; it != config.cluster_steps.end(); ++it) {
			if (ends_with(*it, "_lin") && !ends_with(*(it - 1), "_lin"))
				throw runtime_error("Invalid cluster step sequence: linear steps must be at the start of the list.");
		}
		return config.cluster_steps;
	}
	vector<string> v = { "faster_lin" };
	if (approx_id < 90)
		v.push_back("fast_lin");
	if (approx_id < 40)
		v.push_back("linclust-20_lin");
	else if (approx_id < 80)
		v.push_back("linclust-40_lin");
	if (linear)
		return v;
	if (approx_id < 80)
		v.push_back("default");
	else
		v.push_back("fast");	
	if (approx_id < 50)
		v.push_back("more-sensitive");
	return v;
}

bool is_linclust(const vector<string>& steps) {
	for (const string& s : steps) {
		if (!ends_with(s, "_lin"))
			return false;
	}
	return true;
}

vector<string> default_round_approx_id(int steps) {
	return {};
	/*switch (steps) {
	case 1:
		return {};
	case 2:
		return { "27.0" };
	default:
		return { "27.0", "0.0" };
	}*/
}

vector<string> default_round_cov(int steps) {
	return {};
	/*switch (steps) {
	case 1:
		return {};
	case 2:
		return { "85.0" };
	default:
		return { "87.0", "85.0" };
	}*/
}

static int round_ccd(const string& depth) {
	stringstream ss(depth);
	int i;
	ss >> i;
	if (ss.fail() || !ss.eof())
		throw runtime_error("Invalid number format for --connected-component-depth");
	return i;
}

int round_ccd(std::vector<std::string> param, int round, int round_count, bool linear) {
	if (param.size() > 1 && param.size() != (size_t)round_count)
		throw runtime_error("Parameter count for --connected-component-depth has to be 1 or the number of cascaded clustering rounds.");
	if (param.size() == 0)
		return 0;
	if (param.size() > 1)
		return round_ccd(param[round]);
	return (round == round_count - 1) ^ linear ? 1 : round_ccd(param[0]);
}

}