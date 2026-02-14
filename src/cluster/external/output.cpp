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
#include <sstream>
#include "external.h"
#include "file_array.h"
#include "radix_sort.h"
#include "input_buffer.h"
#include "util/string/string.h"
#include "volume.h"

using std::ostringstream;
using std::runtime_error;
using std::string;
using std::vector;
using std::atomic;
using std::unique_ptr;
using std::ifstream;
using std::pair;
using std::ofstream;

static vector<OId> read_clustering(Job& job, int round) {
	vector<OId> v(job.max_oid + 1);
	vector<OId>::iterator ptr = v.begin();
	for (size_t i = 0; i < job.volumes; ++i) {
		InputFile in(job.base_dir(round) + PATH_SEPARATOR + "clustering" + PATH_SEPARATOR + "volume" + std::to_string(i), InputFile::NO_AUTODETECT);
		const OId n = in.file_size() / sizeof(OId);
		in.read(&*ptr, n);
		ptr += n;
		in.close();
	}
	return v;
}

static vector<OId> merge(Job& job) {
	vector<OId> inner = read_clustering(job, job.round());
	for (int round = job.round() - 1; round >= 0; --round) {
		vector<OId> outer = read_clustering(job, round);
		for (int64_t i = 0; i < (int64_t)outer.size(); ++i)
			outer[i] = inner[outer[i]];
		inner = std::move(outer);
	}
	return inner;
}

static OId output_oids(const vector<OId>& merged) {
	FILE* out = fopen(config.output_file.c_str(), "wt");
	if (!out)
		throw runtime_error("Error opening file: " + config.output_file);
	int64_t n = 0;
	for (int64_t i = 0; i < (int64_t)merged.size(); ++i) {
		if (merged[i] == i)
			++n;
		fprintf(out, "%" PRId64 "\t%" PRId64 "\n", merged[i], i);
	}
	fclose(out);
	return n;
}

struct AccMapping {
	static constexpr bool POD = false;
	static constexpr OId NIL = std::numeric_limits<OId>::max();
	OId rep = NIL, member = NIL;
	std::pmr::string rep_acc, member_acc;
	AccMapping(OId rep, OId member, std::pmr::string&& member_acc, std::pmr::monotonic_buffer_resource& pool) :
		rep(rep),
		member(member),
		rep_acc(&pool),
		member_acc(std::move(member_acc))
	{
	}
	AccMapping(std::pmr::monotonic_buffer_resource& pool):
		rep(NIL),
		member(NIL),
		rep_acc(&pool),
		member_acc(&pool)
	{
	}
	OId key() const {
		return rep;
	}
	bool operator<(const AccMapping& m) const {
		return rep < m.rep || (rep == m.rep && member < m.member);
	}
	friend void serialize(const AccMapping& m, CompressedBuffer& buf) {
		buf.write(m.rep);
		buf.write(m.member);
		buf.write(m.rep_acc.c_str(), m.rep_acc.length() + 1);
		buf.write(m.member_acc.c_str(), m.member_acc.length() + 1);
	}
	friend void deserialize(InputFile& f, AccMapping& m) {
		f.read(&m.rep);
		f.read(&m.member);
		f >> m.rep_acc;
		f >> m.member_acc;
	}
};

static RadixedTable output_accs_round1(Job& job, const vector<OId>& merged, const VolumedFile& volumes) {
	const string base_dir = job.root_dir() + "output" + PATH_SEPARATOR, qpath = base_dir + "queue_round1";
	mkdir(base_dir);
	unique_ptr<FileArray> output_files(new FileArray(base_dir, RADIX_COUNT, job.worker_id(), false));
	Atomic q(qpath);
	SimpleThreadPool pool;

	auto worker = [&](const atomic<bool>& stop) {
		const int shift = std::max(bit_length(job.max_oid) - RADIX_BITS, 0);
		BufferArray buffers(*output_files, RADIX_COUNT);
		std::pmr::monotonic_buffer_resource mem_pool;
		uint64_t v = 0;
		while (!stop.load(std::memory_order_relaxed) && (v = q.fetch_add(), v < job.volumes)) {
			job.log("Building output (round 1). Volume=%lli/%lli", v + 1, job.volumes);
			ifstream acc_file(job.root_dir() + "accessions" + PATH_SEPARATOR + std::to_string(v) + ".txt");
			std::pmr::string line(&mem_pool);
			OId oid = volumes[v].oid_begin;
			while (std::getline(acc_file, line)) {
				AccMapping m(merged[oid], oid, std::move(line), mem_pool);
				buffers.write(m.rep >> shift, m);
				++oid;
			}
		}
		};

	for (int i = 0; i < config.threads_; ++i)
		pool.spawn(worker);
	pool.join_all();
	return output_files->buckets();
}

static OId output_accs(Job& job, const vector<OId>& merged, const VolumedFile& db) {
	const RadixedTable round1 = output_accs_round1(job, merged, db);
	const RadixedTable round1_sorted = radix_sort<AccMapping>(job, round1, std::max(bit_length(job.max_oid) - RADIX_BITS, 0));
	std::pmr::monotonic_buffer_resource pool;
	
	const int digits = ::digits(job.max_oid, 10);
	const std::string base_path = job.root_dir() + "output" + PATH_SEPARATOR, queue_path = base_path + "queue_round2";
	Atomic queue(queue_path);
	int64_t bucket;
	atomic<OId> cluster_count(0);
	while (bucket = queue.fetch_add(), bucket < (int64_t)round1_sorted.size()) {
		VolumedFile file(round1_sorted[bucket]);
		InputBuffer<AccMapping> data(file, pool);		
		job.log("Building output (round 2). Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, round1_sorted.size(), Util::String::format(data.size()).c_str(), Util::String::format(data.byte_size()).c_str());
		if (data.size() == 0)
			continue;
		data.sort();
		const OId oid_begin = data.front().rep, oid_end = data.back().rep + 1;
		const pair<vector<Volume>::const_iterator, vector<Volume>::const_iterator> volumes = db.find(oid_begin, oid_end);
		atomic<int64_t> next(0);
		SimpleThreadPool pool;
		auto worker = [&](const atomic<bool>& stop) {
			OId my_cluster_count = 0;
			int64_t volume;
			auto table_ptr = data.cbegin();			
			while (!stop.load(std::memory_order_relaxed) && (volume = next.fetch_add(1, std::memory_order_relaxed), volume < volumes.second - volumes.first)) {
				unique_ptr<ofstream> output_file;
				const Volume& v = volumes.first[volume];
				while (table_ptr->rep < v.oid_begin)
					++table_ptr;
				const string name = job.root_dir() + "accessions" + PATH_SEPARATOR + std::to_string((volumes.first + volume) - db.begin()) + ".txt";
				ifstream acc_file(name);
				if(!acc_file.good())
					throw runtime_error("Error opening file: " + name);
				OId file_oid = v.oid_begin;
				string acc;
				while (!stop.load(std::memory_order_relaxed) && file_oid < oid_end && std::getline(acc_file, acc)) {
					if (table_ptr->rep > file_oid) {
						++file_oid;
						continue;
					}
					auto begin = table_ptr;
					OId members = 0;
					while (table_ptr != data.cend() && table_ptr->rep == file_oid) {
						if (!output_file) {
							ostringstream name;
							name << config.output_file << '.' << std::setw(digits) << std::setfill('0') << table_ptr->rep;
							output_file.reset(new ofstream(name.str()));
						}
						*output_file << acc << '\t' << table_ptr->member_acc << std::endl; // '\t' << table_ptr->rep << '\t' << table_ptr->member << std::endl;
						++table_ptr;
						++members;
					}
					if (members > 0)
						++my_cluster_count;
					++file_oid;
				}
			}
			cluster_count.fetch_add(my_cluster_count, std::memory_order_relaxed);
			};
		for (int i = 0; i < std::min(config.threads_, int(volumes.second - volumes.first)); ++i)
			pool.spawn(worker);
		pool.join_all();
	}
	return cluster_count.load(std::memory_order_relaxed);
}

void output(Job& job, const VolumedFile& volumes) {
	job.log("Generating output");
	const vector<OId> merged = merge(job);
	const OId n = config.oid_output ? output_oids(merged) : output_accs(job, merged, volumes);
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