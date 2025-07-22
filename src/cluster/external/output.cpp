#include <inttypes.h>
#include "external.h"

using std::vector;
using std::unordered_map;

namespace Cluster {

static vector<int64_t> read_clustering(Job& job, int round) {
	vector<int64_t> v(job.input_count(round));
	vector<int64_t>::iterator it = v.begin();
	for (int i = 0; i < job.volume_count(round); ++i) {
		InputFile in(job.base_dir(round) + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR + "volume" + std::to_string(i));
		const int64_t n = in.file_size() / sizeof(int64_t);
		in.read(&(*it), n);
		it += n;
		in.close();
	}
	return v;
}

void output(Job& job) {
	job.log("Generating output");
	vector<int64_t> inner = read_clustering(job, job.round());
	for (int round = job.round() - 1; round >= 0; --round) {
		VolumedFile reps(job.base_dir(round) + PATH_SEPARATOR + "reps" + PATH_SEPARATOR + "0" + PATH_SEPARATOR + "bucket.tsv");
		vector<int64_t> inner_to_outer(reps.records());
		vector<int64_t>::iterator it = inner_to_outer.begin();
		for (int i = 0; i < (int)reps.size(); ++i) {
			InputFile f(reps[i].path + ".oid");
			f.read(&(*it), reps[i].record_count);
			it += reps[i].record_count;
			f.close();
		}
		unordered_map<int64_t, int64_t> outer_to_inner;
		outer_to_inner.reserve(inner_to_outer.size());
		for (int64_t i = 0; i < (int64_t)inner_to_outer.size(); ++i)
			outer_to_inner[inner_to_outer[i]] = i;
		vector<int64_t> outer = read_clustering(job, round);
		for (int64_t i = 0; i < (int64_t)outer.size(); ++i)
			outer[i] = inner_to_outer[inner[outer_to_inner.at(outer[i])]];
		inner = std::move(outer);
	}
	FILE* out = fopen(config.output_file.c_str(), "wt");
	if (!out)
		throw std::runtime_error("Error opening file: " + config.output_file);
	int64_t n = 0;
	for (int64_t i = 0; i < (int64_t)inner.size(); ++i) {
		if (inner[i] == i)
			++n;
		fprintf(out, PRId64 "\t" PRId64 "\n", inner[i], i);
	}
	fclose(out);
	job.log("Cluster count = %lli", n);
	for (int i = 0; i < job.volume_count(job.round()); ++i)
		remove((job.base_dir(job.round()) + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR + "volume" + std::to_string(i)).c_str());
	for (int round = job.round() - 1; round >= 0; --round) {
		VolumedFile reps(job.base_dir(round) + PATH_SEPARATOR + "reps" + PATH_SEPARATOR + "0" + PATH_SEPARATOR + "bucket.tsv");
		for (int i = 0; i < (int)reps.size(); ++i)
			remove((reps[i].path + ".oid").c_str());
		reps.remove();
		for (int i = 0; i < job.volume_count(round); ++i)
			remove((job.base_dir(round) + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR + "volume" + std::to_string(i)).c_str());
	}
}

}