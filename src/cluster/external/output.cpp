/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

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

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
#include <inttypes.h>
#include "external.h"

using std::vector;

static vector<int64_t> read_clustering(Job& job, int round) {
	vector<int64_t> v(job.max_oid + 1);
	vector<int64_t>::iterator ptr = v.begin();
	for (size_t i = 0; i < job.volumes; ++i) {
		InputFile in(job.base_dir(round) + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR + "volume" + std::to_string(i), InputFile::NO_AUTODETECT);
		const int64_t n = in.file_size() / sizeof(OId);
		in.read(&*ptr, n);
		ptr += n;
		in.close();
	}
	return v;
}

void output(Job& job) {
	job.log("Generating output");
	vector<int64_t> inner = read_clustering(job, job.round());
	for (int round = job.round() - 1; round >= 0; --round) {
		vector<int64_t> outer = read_clustering(job, round);
		for (int64_t i = 0; i < (int64_t)outer.size(); ++i)
			outer[i] = inner[outer[i]];
		inner = std::move(outer);
	}
	FILE* out = fopen(config.output_file.c_str(), "wt");
	if (!out)
		throw std::runtime_error("Error opening file: " + config.output_file);
	int64_t n = 0;
	for (int64_t i = 0; i < (int64_t)inner.size(); ++i) {
		if (inner[i] == i)
			++n;
		fprintf(out, "%" PRId64 "\t%" PRId64 "\n", inner[i], i);
	}
	fclose(out);
	job.log("Cluster count = %lli", n);
	for (size_t i = 0; i < job.volumes; ++i)
		remove((job.base_dir(job.round()) + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR + "volume" + std::to_string(i)).c_str());
	for (int round = job.round() - 1; round >= 0; --round) {
		VolumedFile reps(job.base_dir(round) + PATH_SEPARATOR + "reps" + PATH_SEPARATOR + "reps.tsv");		
		reps.remove();
		for (size_t i = 0; i < job.volumes; ++i)
			remove((job.base_dir(round) + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR + "volume" + std::to_string(i)).c_str());
	}
}