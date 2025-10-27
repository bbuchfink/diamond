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
	if (!config.cluster_steps.empty())
		return config.cluster_steps;
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

int round_ccd(int round, int round_count, bool linear) {
	if (config.connected_component_depth.size() > 1 && config.connected_component_depth.size() != (size_t)round_count)
		throw runtime_error("Parameter count for --connected-component-depth has to be 1 or the number of cascaded clustering rounds.");
	if (config.connected_component_depth.size() == 0)
		return 0;
	if (config.connected_component_depth.size() > 1)
		return round_ccd(config.connected_component_depth[round]);
	return (round == round_count - 1) ^ linear ? 1 : round_ccd(config.connected_component_depth[0]);
}

}