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
#include "util/io/file.h"

using std::vector;
using std::string;
using std::runtime_error;

static vector<OId> read_clusters(const string& path, OId max_oid) {
	File in(path, "rb");
	const size_t count = max_oid + 1;
	if (in.size() != count * sizeof(OId))
		throw runtime_error("Invalid binary clustering file size: " + path);
	vector<OId> v(count);
	in.read(v.data(), v.size() * sizeof(OId));
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
	vector<OId> mapping = read_clusters(job.base_dir(0) + "clusters.bin", job.max_oid());
	rmdir(job.base_dir(0));
	for (int r = 1; r <= job.round(); ++r) {
		chain_round(mapping, job.base_dir(r) + "clusters.bin", job.max_oid());
		rmdir(job.base_dir(r));
	}
	return mapping;
}
