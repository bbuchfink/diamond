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

#include "multinode.h"
#include <inttypes.h>
#include "volume.h"
#include "data/sequence_file.h"
#include "util/parallel/simple_thread_pool.h"

using std::ostringstream;
using std::thread;
using std::vector;
using std::string;
using std::unique_ptr;
using std::atomic;
using std::ofstream;
using std::ifstream;
using std::runtime_error;

string get_reps(Job& job, const VolumedFile& volumes) {
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
			job.log("Writing representatives. Volume=%lli/%lli Records=%s", v + 1, volumes.size(), Util::String::format(volumes[v].record_count).c_str());
			string id_file = job.base_dir() + "rep_ids" + std::to_string(v);
			ifstream rep_ids(id_file);
			if (!rep_ids)
				throw runtime_error("Error opening file " + id_file);
			unique_ptr<SequenceFile> in(SequenceFile::auto_create({ volumes[v].path }));
			string id;
			vector<Letter> seq;
			const string out_file = base_dir + std::to_string(v) + ".faa";
			ofstream out(out_file);
			if (!out)
				throw runtime_error("Error opening file " + out_file);
			OId count = 0, file_oid = volumes[v].oid_begin;
			OId rep;
			rep_ids >> rep;
			while (!stop.load(std::memory_order_relaxed) && in->read_seq(seq, id, nullptr)) {
				if (job.round() > 0)
					file_oid = std::atoll(id.c_str());
				if (rep == file_oid) {
					out << ">"
						<< rep
						<< '\n'
						<< Sequence(seq).to_string()
						<< '\n';
					++count;
					rep_ids >> rep;
				}
				++file_oid;
			}
			if (rep_ids)
				throw runtime_error("Failed to find oid " + std::to_string(rep) + " in file " + volumes[v].path);
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