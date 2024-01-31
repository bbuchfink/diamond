/****
DIAMOND protein aligner
Copyright (C) 2019-2023 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#include <sstream>
#include "cascaded.h"

using std::stringstream;
using std::string;
using std::vector;

namespace Cluster {

vector<string> cluster_steps(double approx_id, bool linear) {
	if (!config.cluster_steps.empty())
		return config.cluster_steps;
	vector<string> v = { "faster_lin" };
	if (linear)
		return v;
	v.push_back("fast");
	if (approx_id < 90)
		v.push_back("default");
	if (approx_id < 50)
		v.push_back("more-sensitive");
	return v;
}

vector<string> default_round_approx_id(int steps) {
	switch (steps) {
	case 1:
		return {};
	case 2:
		return { "27.0" };
	default:
		return { "27.0", "0.0" };
	}
}

vector<string> default_round_cov(int steps) {
	switch (steps) {
	case 1:
		return {};
	case 2:
		return { "85.0" };
	default:
		return { "87.0", "85.0" };
	}
}

static int round_ccd(int round) {
	stringstream ss(config.connected_component_depth[round]);
	int i;
	ss >> i;
	if (ss.fail() || !ss.eof())
		throw std::runtime_error("Invalid number format for --connected-component-depth");
	return i;
}

int round_ccd(int round, int round_count) {
	if (config.connected_component_depth.size() > 1 && config.connected_component_depth.size() != (size_t)round_count)
		throw std::runtime_error("Parameter count for --connected-component-depth has to be 1 or the number of cascaded clustering rounds.");
	if (config.connected_component_depth.size() == 0)
		return 0;
	if (config.connected_component_depth.size() == 1)
		return round_ccd(0);
	return round_ccd(round);
}

}