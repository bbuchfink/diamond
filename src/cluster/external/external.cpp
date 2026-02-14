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
#include <type_traits>
#include <inttypes.h>
#include <cstdarg>
#include <string>
#include "../basic/config.h"
#include "util/tsv/tsv.h"
#include "util/parallel/atomic.h"
#include "util/log_stream.h"
#include "data/sequence_file.h"
#include "basic/reduction.h"
#include "basic/shape_config.h"
#include "search/search.h"
#include "external.h"
#include "radix_sort.h"
#define _REENTRANT
#include "lib/ips4o/ips4o.hpp"
#include "util/algo/hyperloglog.h"
#include "util/sequence//sequence.h"
#include "util/string/string.h"
#include "../cascaded/cascaded.h"
#include "build_pair_table.h"
#include "util/parallel/simple_thread_pool.h"
#include "input_buffer.h"

using std::unique_ptr;
using std::thread;
using std::string;
using std::vector;
using std::endl;
using std::array;
using std::mutex;
using std::lock_guard;
using std::pair;
using std::atomic;
using Util::String::format;
using std::shared_ptr;
using std::to_string;
using std::ofstream;
using std::pmr::monotonic_buffer_resource;

void Job::log(const char* format, ...) {
	char buffer[1024];
	const long long int t = std::chrono::duration_cast<std::chrono::duration<long long int>>(std::chrono::system_clock::now() - start_).count();
	char* ptr = buffer + snprintf(buffer, 1024, "[%" PRId64 ", %lli] ", worker_id_, t);
	va_list args;
	va_start(args, format);
	int i = vsnprintf(ptr, 1024 - (ptr - buffer), format, args);
#ifdef WIN32
	ptr[i++] = '\r';
#endif
	ptr[i++] = '\n';
	ptr[i] = '\0';
	log_stream << buffer;
	log_file_->push(buffer);
	va_end(args);
}

void Job::log(const ClusterStats& stats) {
	std::stringstream ss;
	ss << stats.masking_stat;
	log(ss.str().c_str());
	log("Seeds considered: %" PRIu64, stats.seeds_considered);
	log("Seeds indexed: %" PRIu64, stats.seeds_indexed);	
	log("Extensions computed: %" PRIu64, stats.extensions_computed);
	log("Alignments passing e-value filter: %" PRIu64, stats.hits_evalue_filtered);
	log("Alignments passing all filters: %" PRIu64, stats.hits_filtered);
}

#pragma pack(1)
struct ChunkTableEntry {
	static constexpr bool POD = true;
	ChunkTableEntry():
		oid(),
		chunk()
	{}
	ChunkTableEntry(int64_t oid, int32_t chunk):
		oid(oid),
		chunk(chunk)
	{}
	int64_t key() const {
		return oid;
	}
	bool operator<(const ChunkTableEntry& e) const {
		return oid < e.oid || (oid == e.oid && chunk < e.chunk);
	}
	friend void serialize(const ChunkTableEntry& e, CompressedBuffer& buf) {
		buf.write(e.oid);
		buf.write(e.chunk);
	}
	friend void deserialize(InputFile& in, ChunkTableEntry& e) {
		in.read(&e.oid);
		in.read(&e.chunk);
	}
	OId oid;
	int32_t chunk;
} PACKED_ATTRIBUTE;
#pragma pack()

struct SizeCounter {
	void add(int64_t oid, int32_t len) {
		const int64_t x = oid << 17, n = x + int64_t(len + 63) / 64;
		for (int64_t i = x; i < n; ++i)
			hll.add(i);
	}
	HyperLogLog hll;
};

struct ClusterChunk {
	ClusterChunk(Atomic& next_chunk, const string& chunks_path):
		id((int32_t)next_chunk.fetch_add())
	{
		mkdir(chunks_path + std::to_string(id));
		pairs_out.reset(new OutputFile(chunks_path + std::to_string(id) + PATH_SEPARATOR + "pairs"));
	}
	void write(vector<PairEntryShort>& pairs_buffer, SizeCounter& size) {
		lock_guard<mutex> lock(mtx);
		pairs_out->write(pairs_buffer.size());
		pairs_out->write(pairs_buffer.data(), pairs_buffer.size());
		pairs_buffer.clear();
		this->size.merge(size.hll);
		size.hll = HyperLogLog();
	}
	~ClusterChunk() {
		pairs_out->close();
	}
	const int32_t id;
	unique_ptr<OutputFile> pairs_out;
	HyperLogLog size;
	mutex mtx;
};

static pair<RadixedTable, int> build_chunk_table(Job& job, const RadixedTable& pair_table, int64_t max_oid) {
	const int64_t BUF_SIZE = 4096, shift = std::max(bit_length(max_oid) - RADIX_BITS, 0), max_chunk_size = Util::String::interpret_number(config.linclust_chunk_size) / 64,
		max_processed = std::max(std::min(INT64_C(262144), max_chunk_size / config.threads_ / 16), INT64_C(1)); // ???
	const std::string base_path = job.base_dir() + PATH_SEPARATOR + "chunk_table",
		chunks_path = job.base_dir() + PATH_SEPARATOR + "chunks" + PATH_SEPARATOR;
	mkdir(base_path);
	mkdir(chunks_path);
	unique_ptr<FileArray> output_files(new FileArray(base_path, RADIX_COUNT, job.worker_id(), false));
	Atomic queue(base_path + PATH_SEPARATOR + "queue");
	Atomic next_chunk(base_path + PATH_SEPARATOR + "next_chunk");
	shared_ptr<ClusterChunk> current_chunk(new ClusterChunk(next_chunk, chunks_path));
	int64_t total_pairs = 0, total_distinct_pairs = 0;
	int64_t bucket, buckets_processed = 0;
	mutex mtx;
	while (bucket = queue.fetch_add(), bucket < (int64_t)pair_table.size()) {
		VolumedFile file(pair_table[bucket]);
		InputBuffer<PairEntry> data(file);
		job.log("Building chunk table. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, pair_table.size(), format(data.size()).c_str(), format(data.byte_size()).c_str());
		total_pairs += data.size();
		ips4o::parallel::sort(data.begin(), data.end(), std::less<PairEntry>(), config.threads_);
		SimpleThreadPool pool;
		auto worker = [&](const atomic<bool>& stop, int thread_id) {
			shared_ptr<ClusterChunk> my_chunk(current_chunk);
			BufferArray buffers(*output_files, RADIX_COUNT);
			vector<PairEntryShort> pairs_buffer;
			SizeCounter size;
			int64_t distinct_pairs = 0, processed = 0;
			auto it = merge_keys(data.begin(thread_id), data.end(thread_id), PairEntry::Key());
			while (!stop.load(std::memory_order_relaxed) && it.good()) {
				const uint64_t rep_oid = it.begin()->rep_oid;
				buffers.write(rep_oid >> shift, ChunkTableEntry(rep_oid, my_chunk->id));
				size.add(rep_oid, it.begin()->rep_len);
				processed += it.begin()->rep_len;
				for (auto j = it.begin(); j < it.end(); ++j) {
					if (j > it.begin() && j->member_oid == (j - 1)->member_oid)
						continue;
					const int radix = int(j->member_oid >> shift);
					buffers.write(radix, ChunkTableEntry(j->member_oid, my_chunk->id));
					size.add(j->member_oid, j->member_len);
					pairs_buffer.emplace_back(rep_oid, j->member_oid);
					++distinct_pairs;
					processed += j->member_len;
					if (processed >= max_processed) {
						my_chunk->write(pairs_buffer, size);
						processed = 0;
						bool new_chunk = false;
						if (my_chunk != current_chunk) {
							my_chunk = current_chunk;
							new_chunk = true;
						}
						else {
							const int64_t est = my_chunk->size.estimate();
							if (est >= max_chunk_size) {
								lock_guard<mutex> lock(mtx);
								if (my_chunk == current_chunk) {
									log_stream << "build_chunk_table chunk=" << current_chunk->id << " est_size=" << est * 64 << endl;
									current_chunk.reset(new ClusterChunk(next_chunk, chunks_path));
									my_chunk = current_chunk;
									new_chunk = true;
								}
								else {
									// should not happen?
								}
							}
						}
						if (new_chunk) {
							buffers.write(int(rep_oid >> shift), ChunkTableEntry(rep_oid, my_chunk->id));
							size.add(rep_oid, it.begin()->rep_len);
							processed += it.begin()->rep_len;
						}
					}
				}
				++it;
			}
			my_chunk->write(pairs_buffer, size);
			total_distinct_pairs += distinct_pairs;
			};
		for (int i = 0; i < data.parts(); ++i)
			pool.spawn(worker, i);
		pool.join_all();
		const int64_t est = current_chunk->size.estimate();
		if (est >= max_chunk_size) {
			log_stream << "build_chunk_table chunk=" << current_chunk->id << " est_size=" << est * 64 << endl;
			current_chunk.reset(new ClusterChunk(next_chunk, chunks_path));
		}
		file.remove();
		++buckets_processed;
	}
	log_stream << "build_chunk_table chunk=" << current_chunk->id << " est_size=" << current_chunk->size.estimate() << " total_pairs=" << total_pairs << " total_distinct_pairs=" << total_distinct_pairs << endl;
	const RadixedTable buckets = output_files->buckets();
	TaskTimer timer("Closing the output files");
	output_files.reset();
	current_chunk.reset();
	timer.go("Waiting for other workers");
	Atomic finished(base_path + PATH_SEPARATOR + "finished");
	finished.fetch_add(buckets_processed);
	finished.await(pair_table.size());
	return { buckets, (int)next_chunk.get() };
}

static void build_chunks(Job& job, const VolumedFile& db, const RadixedTable& chunk_table, int chunk_count) {
	const int64_t BUF_SIZE = 64 * 1024;
	const std::string base_path = job.base_dir() + PATH_SEPARATOR + "chunks" + PATH_SEPARATOR,
		queue_path = base_path + "queue";
	unique_ptr<FileArray> output_files(new FileArray(base_path, chunk_count, job.worker_id(), false, 1024 * 1024 * 1024));
	Atomic queue(queue_path);
	int64_t bucket, buckets_processed = 0;
	atomic<int64_t> oid_counter(0), distinct_oid_counter(0);
	while (bucket = queue.fetch_add(), bucket < (int64_t)chunk_table.size()) {
		VolumedFile file(chunk_table[bucket]);
		InputBuffer<ChunkTableEntry> data(file);
		job.log("Building chunks. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, chunk_table.size(), format(data.size()).c_str(), format(data.byte_size()).c_str());
		if (data.size() > 0) {
			ips4o::parallel::sort(data.begin(), data.end(), std::less<ChunkTableEntry>(), config.threads_);
			const OId oid_begin = data.front().oid, oid_end = data.back().oid + 1;
			const pair<vector<Volume>::const_iterator, vector<Volume>::const_iterator> volumes = db.find(oid_begin, oid_end);
			atomic<int64_t> next(0);
			SimpleThreadPool pool;
			auto worker = [&](const atomic<bool>& stop) {
				int64_t volume;
				auto table_ptr = data.cbegin();
				BufferArray output_bufs(*output_files, chunk_count);
				TextBuffer buf;
				while (!stop.load(std::memory_order_relaxed) && (volume = next.fetch_add(1, std::memory_order_relaxed), volume < volumes.second - volumes.first)) {
					const Volume& v = volumes.first[volume];
					while (table_ptr->oid < v.oid_begin)
						++table_ptr;
					unique_ptr<SequenceFile> in(SequenceFile::auto_create({ v.path }));
					string id;
					vector<Letter> seq;
					OId file_oid = v.oid_begin;
					while (!stop.load(std::memory_order_relaxed) && file_oid < oid_end && in->read_seq(seq, id, nullptr)) {
						if(job.round() > 0)
							file_oid = atoll(id.c_str());
						if (table_ptr->oid > file_oid) {
							++file_oid;
							continue;
						}
						Util::Seq::format(seq, std::to_string(file_oid).c_str(), nullptr, buf, "fasta", amino_acid_traits);
						auto begin = table_ptr;
						while (table_ptr < data.end() && table_ptr->oid == file_oid) {
							if (table_ptr == begin || table_ptr->chunk != table_ptr[-1].chunk) {
								output_bufs.write(table_ptr->chunk, buf.data(), buf.size());
								oid_counter.fetch_add(1, std::memory_order_relaxed);
							}
							++table_ptr;
						}
						buf.clear();
						distinct_oid_counter.fetch_add(1, std::memory_order_relaxed);
						++file_oid;
					}
					in->close();
				}
				};
			for (int i = 0; i < std::min(config.threads_, int(volumes.second - volumes.first)); ++i)
				pool.spawn(worker);
			pool.join_all();
		}
		file.remove();
		++buckets_processed;
	}
	TaskTimer timer("Closing the output files");
	output_files.reset();
	timer.go("Waiting for other workers");
	Atomic finished(base_path + "finished");
	finished.fetch_add(buckets_processed);
	finished.await(chunk_table.size());
	timer.finish();
	log_stream << "build_chunks oids=" << oid_counter << '/' << db.sparse_records() << " distinct_oids=" << distinct_oid_counter << endl;
}

string round(Job& job, const VolumedFile& volumes) {
	::shapes = ShapeConfig(Search::shape_codes.at(config.sensitivity), 0);
	if (config.mutual_cover.present()) {
		config.min_length_ratio = config.sensitivity < Sensitivity::LINCLUST_40 ?
			std::min(config.mutual_cover.get_present() / 100 + 0.05, 1.0)
			: config.mutual_cover.get_present() / 100 - 0.05;
	}
	job.log("Starting round %i sensitivity %s %i shapes\n", job.round(), to_string(config.sensitivity).c_str(), ::shapes.count());
	job.set_round(volumes.sparse_records());
	const string pair_table_base = job.base_dir() + PATH_SEPARATOR + "pair_table";
	mkdir(pair_table_base);
	unique_ptr<FileArray> pair_table_files(new FileArray(pair_table_base, RADIX_COUNT, job.worker_id(), false));
	RadixedTable pair_table;
	for (int shape = 0; shape < ::shapes.count(); ++shape) {
		const RadixedTable buckets = build_seed_table(job, volumes, shape);
		const RadixedTable sorted_seed_table = radix_sort<SeedEntry>(job, buckets, 64 - RADIX_BITS);
		pair_table = build_pair_table(job, sorted_seed_table, shape, volumes.max_oid(), *pair_table_files);
	}
	pair_table_files.reset();
	const RadixedTable sorted_pair_table = radix_sort<PairEntry>(job, pair_table, 64 - RADIX_BITS);
	const pair<RadixedTable, int> chunk_table = build_chunk_table(job, sorted_pair_table, volumes.max_oid());
	const RadixedTable sorted_chunk_table = radix_sort<ChunkTableEntry>(job, chunk_table.first, bit_length(volumes.max_oid()) - RADIX_BITS);
	build_chunks(job, volumes, sorted_chunk_table, chunk_table.second);
	const RadixedTable edges = align(job, chunk_table.second, volumes.max_oid());
	if (config.mutual_cover.present()) {
		return cluster_bidirectional(job, edges, volumes);
	}
	else {
		const RadixedTable sorted_edges = radix_sort<Edge>(job, edges, 64 - RADIX_BITS);
		//const vector<string> sorted_edges = read_list(job.base_dir() + PATH_SEPARATOR + "alignments" + PATH_SEPARATOR + "radix_sort_out");
		return cluster(job, sorted_edges, volumes);
	}
}

void external() {
	if (config.output_file.empty())
		throw std::runtime_error("Option missing: output file (--out/-o)");
	config.file_buffer_size = 64 * 1024;
	TaskTimer total;
	VolumedFile volumes(config.database.get_present());
	Job job(volumes.max_oid(), volumes.size());
	if (job.worker_id() == 0) {
		if(config.mutual_cover.present())
			job.log("Bi-directional coverage = %f", config.mutual_cover.get_present());
		else
			job.log("Uni-directional coverage = %f", config.member_cover.get(80));
		job.log("Approx. id = %f", config.approx_min_id.get(0));
		job.log("#Volumes = %lli", volumes.size());
		job.log("#Sequences = %lli", volumes.sparse_records());
	}
	if (config.mutual_cover.present()) {
		config.query_or_target_cover = 0;
		config.query_cover = config.mutual_cover.get_present();
		config.subject_cover = config.mutual_cover.get_present();
	}
	else {
		config.query_or_target_cover = config.member_cover.get(80);
		config.query_cover = 0;
		config.subject_cover = 0;
	}
#ifdef WIN32
	_setmaxstdio(8192);
#endif
	vector<string> steps = Cluster::cluster_steps(config.approx_min_id.get(0), true);
	string reps;
	job.set_round_count((int)steps.size());
	for (size_t i = 0; i < steps.size(); ++i) {
		config.sensitivity = from_string<Sensitivity>(rstrip(steps[i], "_lin"));
		reps = round(job, i == 0 ? volumes : VolumedFile(reps));
		if (i < steps.size() - 1)
			job.next_round();
	}
	Atomic output_lock(job.root_dir() + PATH_SEPARATOR + "output_lock");
	if(output_lock.fetch_add() == 0)
		output(job, volumes);
	log_stream << "Total time = " << (double)total.milliseconds() / 1000 << 's' << endl;
}