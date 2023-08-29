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

#include "cascaded.h"

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

}