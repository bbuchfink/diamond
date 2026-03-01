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

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
#include <numeric>
#include <inttypes.h>
#include <sstream>
#include "basic/config.h"
#include "external.h"
#include "util/parallel/atomic.h"
#define _REENTRANT
#include "lib/ips4o/ips4o.hpp"
#include "util/util.h"
#include "data/sequence_file.h"
#include "util/parallel/simple_thread_pool.h"
#include "file_array.h"
#include "input_buffer.h"

using std::vector;
using std::string;
using std::unique_ptr;
using std::thread;
using std::endl;
using Util::String::format;
using std::atomic;
using std::ofstream;
using std::ostringstream;

static void compute_closure(Job& job, const VolumedFile& volumes, vector<uint64_t>& rep) {
	Partition<int64_t> parts(volumes.max_oid() + 1, config.threads_);
	auto closure_worker = [&](int thread_id) {
		int64_t clusters = 0;
		for (int64_t i = parts.begin(thread_id); i < parts.end(thread_id); ++i) {
			int64_t r = rep[i];
			if (r == i)
				++clusters;
			while (rep[r] != r)
				r = rep[r];
			if (r != rep[i])
				rep[i] = r;
		}
	};
	vector<thread> workers;
	for (int i = 0; i < std::min(config.threads_, (int)parts.parts); ++i)
		workers.emplace_back(closure_worker, i);
	for (auto& t : workers)
		t.join();
	const string output_dir = job.base_dir() + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR;
	mkdir(output_dir);
	for (size_t v = 0; v < volumes.size(); ++v) {
		OutputFile out(output_dir + "volume" + std::to_string(v));
		out.write(&rep[volumes[v].oid_begin], volumes[v].oid_range());
		out.close();
	}
}

static void compute_closure(Job& job, const string& assignment_file, const VolumedFile& volumes) {
	job.log("Computing transitive closure");
	TaskTimer timer("Getting assignment volumes");
	VolumedFile v(assignment_file);
	timer.go("Initializing mapping vector");
	vector<OId> rep(volumes.max_oid() + 1);
	std::iota(rep.begin(), rep.end(), 0);
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
	compute_closure(job, volumes, rep);
	v.remove();
}

static string get_reps(Job& job, const VolumedFile& volumes) {
	if (job.last_round())
		return string();
	const string base_dir = job.base_dir() + PATH_SEPARATOR + "reps" + PATH_SEPARATOR, qpath = base_dir + "queue";
	const string clustering_dir = job.base_dir() + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR;
	mkdir(base_dir);
	unique_ptr<FileStack> reps_list(new FileStack(base_dir + "reps.tsv"));
	
	Atomic q(qpath);
	atomic<int> volumes_processed(0);
	atomic<OId> cluster_count(0);
	vector<thread> workers;
	auto worker = [&](const atomic<bool>& stop, int thread_id) {
		int64_t v = 0;
		vector<Letter> buf;
		while (v = q.fetch_add(), !stop.load(std::memory_order_relaxed) && v < (int64_t)volumes.size()) {
			job.log("Writing representatives. Volume=%lli/%lli Records=%s", v + 1, volumes.size(), format(volumes[v].record_count).c_str());
			fflush(stdout);
			vector<OId> rep(volumes[v].oid_range());
			InputFile clustering_file(clustering_dir + "volume" + std::to_string(v), InputFile::NO_AUTODETECT);
			clustering_file.read(rep.data(), rep.size());
			clustering_file.close();
			unique_ptr<SequenceFile> in(SequenceFile::auto_create({ volumes[v].path }));
			string id;
			vector<Letter> seq;
			OId file_oid = volumes[v].oid_begin, table_oid = file_oid;
			vector<OId>::const_iterator rep_it = rep.begin();
			const string out_file = base_dir + std::to_string(v) + ".faa";
			ofstream out(out_file);
			OId count = 0, max_oid = 0, min_oid = std::numeric_limits<OId>::max();
			while (!stop.load(std::memory_order_relaxed) && in->read_seq(seq, id, nullptr)) {
				if (job.round() > 0) {
					file_oid = atoll(id.c_str());
				}
				while(table_oid < file_oid) {
					++table_oid;
					++rep_it;
				}
				if (*rep_it == file_oid) {
					out << ">"
						<< file_oid
						<< '\n'
						<< Sequence(seq).to_string()
						<< '\n';
					++count;
					max_oid = std::max(max_oid, file_oid);
					min_oid = std::min(min_oid, file_oid);
				}
				++file_oid;
				++rep_it;
				++table_oid;
			}
			in->close();
			out.close();
			ostringstream ss;
			ss << out_file << '\t' << count << '\t' << volumes[v].oid_begin << '\t' << volumes[v].oid_end << '\n';
			reps_list->push(ss.str());
			volumes_processed.fetch_add(1, std::memory_order_relaxed);
			cluster_count.fetch_add(count, std::memory_order_relaxed);
		}
		};
	SimpleThreadPool pool;
	for (int i = 0; i < config.threads_; ++i)
		pool.spawn(worker, i);
	int n = 0;
	pool.join_all();
	job.log("Representatives written: %" PRIu64, cluster_count.load());
	TaskTimer timer("Closing the output files");
	Atomic finished(base_dir + "finished");
	finished.fetch_add(volumes_processed);
	finished.await((int)volumes.size());
	const string out = reps_list->file_name();
	reps_list.reset();
	return out;
}

string cluster(Job& job, const RadixedTable& edges, const VolumedFile& volumes) {
	const std::string base_path = job.base_dir() + PATH_SEPARATOR + "alignments",
		queue_path = base_path + PATH_SEPARATOR + "queue",
		clustering_path = job.base_dir() + PATH_SEPARATOR + "clustering";
	mkdir(clustering_path);
	unique_ptr<FileArray> output_file(new FileArray(clustering_path, 1, job.worker_id(), false));
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
		file.remove();
		++buckets_processed;
	}
	TaskTimer timer("Closing the output files");
	const string assignment_file = output_file->bucket(0);
	output_file.reset();
	Atomic finished(clustering_path + PATH_SEPARATOR + "finished");
	const int64_t f = finished.fetch_add(buckets_processed);
	Atomic closure_finished(clustering_path + PATH_SEPARATOR + "closure_finished");
	if (f + buckets_processed < (int)edges.size())
		closure_finished.await(1);
	else {
		compute_closure(job, assignment_file, volumes);
		closure_finished.fetch_add();
	}
	return get_reps(job, volumes);
}

string cluster_bidirectional(Job& job, const RadixedTable& edges, const VolumedFile& volumes) {
	Atomic lock(job.base_dir() + PATH_SEPARATOR + "cluster_bidirectional_lock");
	Atomic finished(job.base_dir() + PATH_SEPARATOR + "cluster_bidirectional_finished");
	if (lock.fetch_add() == 0) {
		job.log("Computing clustering (bi-directional coverage)");
		vector<uint32_t> degree(volumes.max_oid() + 1, 0);
		for (int bucket = 0; bucket < (int)edges.size(); ++bucket) {
			VolumedFile file(edges[bucket]);
			InputBuffer<Edge> data(file);
			job.log("Getting node degrees. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, edges.size(), format(data.size()).c_str(), format(data.byte_size()).c_str());
			const auto end = data.cend();
			for (auto i = data.cbegin(); i < end; ++i) {
				++degree[i->member_oid];
				++degree[i->rep_oid];
			}
		}
		vector<uint64_t> rep(volumes.max_oid() + 1);
		std::iota(rep.begin(), rep.end(), (int64_t)0);
		for (int bucket = 0; bucket < (int)edges.size(); ++bucket) {
			VolumedFile file(edges[bucket]);
			InputBuffer<Edge> data(file);
			job.log("Assigning reps. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, edges.size(), format(data.size()).c_str(), format(data.byte_size()).c_str());
			const auto end = data.cend();
			for (auto i = data.cbegin(); i < end; ++i) {
				if (degree[i->rep_oid] > degree[rep[i->member_oid]] || (degree[i->rep_oid] == degree[rep[i->member_oid]] && i->rep_oid < rep[i->member_oid]))
					rep[i->member_oid] = i->rep_oid;
				if (degree[i->member_oid] > degree[rep[i->rep_oid]] || (degree[i->member_oid] == degree[rep[i->rep_oid]] && i->member_oid < rep[i->rep_oid]))
					rep[i->rep_oid] = i->member_oid;
			}
		}
		compute_closure(job, volumes, rep);
		finished.fetch_add();
		for (int bucket = 0; bucket < (int)edges.size(); ++bucket)
			VolumedFile(edges[bucket]).remove();
	}
	else
		finished.await(1);
	return get_reps(job, volumes);
}