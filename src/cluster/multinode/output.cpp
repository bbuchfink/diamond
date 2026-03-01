/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#include <inttypes.h>
#include <fstream>
#include <limits>
#include "multinode.h"
#include "volume.h"
#include "basic/config.h"

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
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
	return v;
}

static void chain_round(vector<OId>& mapping, const string& path, OId max_oid) {
	const vector<OId> next = read_clusters(path, max_oid);
	for (OId& c : mapping)
		if (c != std::numeric_limits<OId>::max())
			c = next[c];
}

static vector<OId> build_merged(Job& job) {
	vector<OId> mapping = read_clusters(job.base_dir(0) + PATH_SEPARATOR + "clusters.tsv", job.max_oid);	
	for (int r = 1; r <= job.round(); ++r)
		chain_round(mapping, job.base_dir(r) + PATH_SEPARATOR + "clusters.tsv", job.max_oid);
	return mapping;
}

static OId output_oids(Job& job, const vector<OId>& merged) {
	FILE* out = fopen(config.output_file.c_str(), "wt");
	if (!out)
		throw runtime_error("Error opening output file: " + config.output_file);
	OId n = 0;
	for (OId i = 0; i <= job.max_oid; ++i) {
		if (merged[i] == i) ++n;
		fprintf(out, "%" PRId64 "\t%" PRId64 "\n", (int64_t)merged[i], (int64_t)i);
	}
	fclose(out);
	return n;
}

static OId output_accs(Job& job, const vector<OId>& merged, const VolumedFile& volumes) {
	vector<string> acc(job.max_oid + 1);
	for (size_t v = 0; v < volumes.size(); ++v) {
		const string path = job.root_dir() + "accessions" + PATH_SEPARATOR + std::to_string(v) + ".txt";
		ifstream in(path);
		if (!in.good())
			throw runtime_error("Error opening accessions file: " + path);
		OId oid = volumes[v].oid_begin;
		while (oid < volumes[v].oid_end && std::getline(in, acc[oid]))
			++oid;
	}
	ofstream out(config.output_file);
	if (!out.good())
		throw runtime_error("Error opening output file: " + config.output_file);
	OId n = 0;
	for (OId i = 0; i <= job.max_oid; ++i) {
		const OId centroid = merged[i];
		out << acc[centroid] << '\t' << acc[i] << '\n';
		if (centroid == i) ++n;
	}
	return n;
}

void merge(Job& job, const VolumedFile& volumes) {
	job.log("Merging clusterings");
	const vector<OId> merged = build_merged(job);
	const OId n = config.oid_output
		? output_oids(job, merged)
		: output_accs(job, merged, volumes);
	job.log("Total clusters: %" PRId64, (int64_t)n);
}