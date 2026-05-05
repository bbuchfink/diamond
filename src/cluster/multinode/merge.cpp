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

#include "multinode.h"

using std::vector;
using std::ifstream;
using std::string;
using std::runtime_error;

static vector<OId> read_clusters(const string& path, OId max_oid) {
	constexpr OId NIL = std::numeric_limits<OId>::max();
	vector<OId> v(max_oid + 1, NIL);
	ifstream in(path);
	if (!in.good())
		throw runtime_error("Error opening clustering file: " + path);
	OId rep, member;
	while (in >> rep >> member)
		v[member] = rep;
	in.close();
	remove_tmp_file(path);
	return v;
}

static void chain_round(vector<OId>& mapping, const string& path, OId max_oid) {
	const vector<OId> next = read_clusters(path, max_oid);
	for (OId& c : mapping)
		if (c != std::numeric_limits<OId>::max())
			c = next[c];
}

vector<OId> build_merged(Job& job) {
	vector<OId> mapping = read_clusters(job.base_dir(0) + PATH_SEPARATOR + "clusters.tsv", job.max_oid());
	rmdir(job.base_dir(0));
	for (int r = 1; r <= job.round(); ++r) {
		chain_round(mapping, job.base_dir(r) + PATH_SEPARATOR + "clusters.tsv", job.max_oid());
		rmdir(job.base_dir(r));
	}
	return mapping;
}