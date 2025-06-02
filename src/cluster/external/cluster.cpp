#include <numeric>
#include "basic/config.h"
#include "external.h"
#include "util/parallel/atomic.h"
#define _REENTRANT
#include "lib/ips4o/ips4o.hpp"
#include "util/util.h"

using std::vector;
using std::string;
using std::unique_ptr;
using std::thread;
using std::endl;
using Util::String::format;
using std::atomic;

namespace Cluster {

struct Assignment {
	int64_t member_oid, rep_oid;
};

void cluster(Job& job, const vector<string>& edges, int64_t db_size) {
	const std::string base_path = job.base_dir() + PATH_SEPARATOR + "alignments",
		queue_path = base_path + PATH_SEPARATOR + "queue",
		clustering_path = config.tmpdir + PATH_SEPARATOR + "clustering";
	mkdir(clustering_path);
	unique_ptr<FileArray> output_file(new FileArray(clustering_path, 1, job.worker_id()));
	Atomic queue(queue_path);
	int64_t bucket, buckets_processed = 0;
	while (bucket = queue.fetch_add(), bucket < (int64_t)edges.size()) {
		VolumedFile file(edges[bucket]);
		InputBuffer<Edge> data(file);
		job.log("Clustering. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, edges.size(), format(data.size()).c_str(), format(data.byte_size()).c_str());
		ips4o::parallel::sort(data.begin(), data.end(), std::less<Edge>(), config.threads_);
		auto worker = [&](int thread_id) {
			BufferArray buffer(*output_file, 1);
			auto it = merge_keys(data.begin(thread_id), data.end(thread_id), Edge::Member());
			while (it.good()) {
				if (it.begin()->member_len < it.begin()->rep_len || (it.begin()->member_len == it.begin()->rep_len && it.begin()->member_oid > it.begin()->rep_oid))
					buffer.write(0, Assignment{ it.begin()->member_oid, it.begin()->rep_oid });
				++it;
			}
			};
		vector<thread> workers;
		for (int i = 0; i < data.parts(); ++i)
			workers.emplace_back(worker, i);
		for (auto& t : workers)
			t.join();
		++buckets_processed;
	}
	TaskTimer timer("Closing the output files");
	const string assignment_file = output_file->bucket(0);
	output_file.reset();
	Atomic finished(clustering_path + PATH_SEPARATOR + "finished");
	const int f = finished.fetch_add(buckets_processed);
	if (f + buckets_processed < (int)edges.size())
		return;
	job.log("Computing transitive closure");
	timer.go("Getting assignment volumes");
	VolumedFile v(assignment_file);
	timer.go("Initializing mapping vector");
	vector<int64_t> rep(db_size);
	std::iota(rep.begin(), rep.end(), (int64_t)0);
	timer.go("Reading assignments");
	atomic<int> next(0);
	auto read_worker = [&](int thread_id) {
		int i;
		while (i = next.fetch_add(1, std::memory_order_relaxed), i < (int)v.size()) {
			job.log("Reading assignments thread_id=%i volume=%i", thread_id, i);
			InputFile f(v[i].path);
			Assignment a;
			while (f.read(&a, 1) == 1) {
				rep[a.member_oid] = a.rep_oid;
			}
			f.close();
		}
		};
	vector<thread> workers;
	for (int i = 0; i < std::min(config.threads_, (int)v.size()); ++i)
		workers.emplace_back(read_worker, i);
	for (auto& t : workers)
		t.join();
	workers.clear();
	timer.go("Computing closure");
	Partition<int64_t> parts(db_size, config.threads_);
	auto closure_worker = [&](int thread_id) {
		for (int64_t i = parts.begin(thread_id); i < parts.end(thread_id); ++i) {
			int64_t r = rep[i];
			while (rep[r] != r)
				r = rep[r];
			if (r != rep[i])
				rep[i] = r;
		}
		};
	for (int i = 0; i < std::min(config.threads_, (int)parts.parts); ++i)
		workers.emplace_back(closure_worker, i);
	for (auto& t : workers)
		t.join();
	timer.go("Generating output");
	FILE* out = fopen(config.output_file.c_str(), "wt");
	if (!out)
		throw std::runtime_error("Error opening file: " + config.output_file);
	int64_t n = 0;
	for (int64_t i = 0; i < db_size; ++i) {
		if (rep[i] == i)
			++n;
#ifdef _MSC_VER
		fprintf(out, "%lli\t%lli\n", rep[i], i);
#else
		fprintf(out, "%li\t%li\n", rep[i], i);
#endif
	}
	timer.finish();
	fclose(out);
	job.log("Cluster count = %l", n);
}

}