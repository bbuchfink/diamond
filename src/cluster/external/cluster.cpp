/****
Copyright © 2013-2025 Benjamin J. Buchfink

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

static void compute_closure(Job& job, const VolumedFile& volumes, vector<int64_t>& rep) {
	Partition<int64_t> parts(volumes.records(), config.threads_);
	atomic<int64_t> cluster_count(0);
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
		cluster_count.fetch_add(clusters, std::memory_order_relaxed);
		};
	vector<thread> workers;
	for (int i = 0; i < std::min(config.threads_, (int)parts.parts); ++i)
		workers.emplace_back(closure_worker, i);
	for (auto& t : workers)
		t.join();
	job.log("Cluster count = %lli", (int64_t)cluster_count);
	const string output_dir = job.base_dir() + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR;
	mkdir(output_dir);
	for (size_t v = 0; v < volumes.size(); ++v) {
		OutputFile out(output_dir + "volume" + std::to_string(v));
		out.write(&rep[volumes[v].oid_begin], volumes[v].record_count);
		out.close();
	}
}

static void compute_closure(Job& job, const string& assignment_file, const VolumedFile& volumes) {
	job.log("Computing transitive closure");
	TaskTimer timer("Getting assignment volumes");
	VolumedFile v(assignment_file);
	timer.go("Initializing mapping vector");
	vector<int64_t> rep(volumes.records());
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
	compute_closure(job, volumes, rep);
	v.remove();
}

static string get_reps(Job& job, const VolumedFile& volumes) {
	if (job.last_round())
		return string();
	const string base_dir = job.base_dir() + PATH_SEPARATOR + "reps" + PATH_SEPARATOR, qpath = base_dir + "queue";
	const string clustering_dir = job.base_dir() + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR;
	mkdir(base_dir);
	unique_ptr<FileArray> output_files(new FileArray(base_dir, 1, job.worker_id()));
	
	Atomic q(qpath);
	atomic<int> volumes_processed(0);
	vector<thread> workers;
	auto worker = [&](int thread_id) {
		BufferArray buffers(*output_files, 1);
		int64_t v = 0;
		vector<Letter> buf;
		while (v = q.fetch_add(), v < (int64_t)volumes.size()) {
			job.log("Writing representatives. Volume=%lli/%lli Records=%s", v + 1, volumes.size(), format(volumes[v].record_count).c_str());
			vector<int64_t> rep(volumes[v].record_count);
			InputFile clustering_file(clustering_dir + "volume" + std::to_string(v));
			clustering_file.read(rep.data(), volumes[v].record_count);
			clustering_file.close();
			unique_ptr<SequenceFile> in(SequenceFile::auto_create({ volumes[v].path }));
			string id;
			vector<Letter> seq;
			int64_t oid = volumes[v].oid_begin;
			vector<int64_t>::const_iterator rep_it = rep.begin();
			string fasta_record;
			while (in->read_seq(seq, id, nullptr)) {
				if (*rep_it == oid) {
					fasta_record = ">";
					fasta_record += std::to_string(oid);
					fasta_record += '\n';
					fasta_record += Sequence(seq).to_string();
					fasta_record += '\n';
					buffers.write(0, fasta_record.data(), fasta_record.length(), 1);
				}
				++oid;
				++rep_it;
			}
			in->close();
			volumes_processed.fetch_add(1, std::memory_order_relaxed);
		}
		};
	for (int i = 0; i < config.threads_; ++i)
		workers.emplace_back(worker, i);
	for (auto& t : workers)
		t.join();
	const vector<string> buckets = output_files->buckets();
	TaskTimer timer("Closing the output files");
	output_files.reset();
	Atomic finished(base_dir + "finished");
	finished.fetch_add(volumes_processed);
	finished.await((int)volumes.size());
	return buckets.front();
}

string cluster(Job& job, const vector<string>& edges, const VolumedFile& volumes) {
	const int64_t db_size = volumes.records();
	const std::string base_path = job.base_dir() + PATH_SEPARATOR + "alignments",
		queue_path = base_path + PATH_SEPARATOR + "queue",
		clustering_path = job.base_dir() + PATH_SEPARATOR + "clustering";
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
		file.remove();
		++buckets_processed;
	}
	TaskTimer timer("Closing the output files");
	const string assignment_file = output_file->bucket(0);
	output_file.reset();
	Atomic finished(clustering_path + PATH_SEPARATOR + "finished");
	const int f = finished.fetch_add(buckets_processed);
	Atomic closure_finished(clustering_path + PATH_SEPARATOR + "closure_finished");
	if (f + buckets_processed < (int)edges.size())
		closure_finished.await(1);
	else {
		compute_closure(job, assignment_file, volumes);
		closure_finished.fetch_add();
	}
	return get_reps(job, volumes);
}

string cluster_bidirectional(Job& job, const vector<string>& edges, const VolumedFile& volumes) {
	Atomic lock(job.base_dir() + PATH_SEPARATOR + "cluster_bidirectional_lock");
	Atomic finished(job.base_dir() + PATH_SEPARATOR + "cluster_bidirectional_finished");
	if (lock.fetch_add() == 0) {
		job.log("Computing clustering (bi-directional coverage)");
		vector<uint32_t> degree(volumes.records(), 0);
		for (int bucket = 0; bucket < (int)edges.size(); ++bucket) {
			VolumedFile file(edges[bucket]);
			InputBuffer<Edge> data(file);
			job.log("Getting node degrees. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, edges.size(), format(data.size()).c_str(), format(data.byte_size()).c_str());
			const Edge* end = data.end();
			for (const Edge* i = data.begin(); i < end; ++i) {
				++degree[i->member_oid];
				++degree[i->rep_oid];
			}
		}
		vector<int64_t> rep(volumes.records());
		std::iota(rep.begin(), rep.end(), (int64_t)0);
		for (int bucket = 0; bucket < (int)edges.size(); ++bucket) {
			VolumedFile file(edges[bucket]);
			InputBuffer<Edge> data(file);
			job.log("Assigning reps. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, edges.size(), format(data.size()).c_str(), format(data.byte_size()).c_str());
			const Edge* end = data.end();
			for (const Edge* i = data.begin(); i < end; ++i) {
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

}