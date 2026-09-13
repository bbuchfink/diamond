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

#include <inttypes.h>
#include <fstream>
#include <unordered_map>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <climits>
#include "multinode.h"
#include "volume.h"
#include "basic/config.h"
#include "util/memory/memory_resource.h"
#include "util/io/compressed_buffer.h"
#include "file_array.h"
#include "input_buffer.h"
#include "cluster.h"
#include "util/parallel/simple_thread_pool.h"
#include "data/sequence_file.h"
#include "util/sequence/sequence.h"
#include "util/text_buffer.h"

using std::endl;
using std::unordered_map;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::runtime_error;
using std::unique_ptr;
using std::atomic;

static OId output_oids(Job& job, const vector<OId>& merged) {
	FILE* out = fopen(config.output_file.c_str(), "wt");
	if (!out)
		throw runtime_error("Error opening output file: " + config.output_file);
	OId n = 0;
	for (OId i = 0; i <= job.max_oid(); ++i) {
		if (merged[i] == i) ++n;
		fprintf(out, "%" PRId64 "\t%" PRId64 "\n", (int64_t)merged[i], (int64_t)i);
	}
	fclose(out);
	return n;
}

struct AccMapping {
	static constexpr bool POD = false;
	static constexpr OId NIL = std::numeric_limits<OId>::max();
	OId rep = NIL, member = NIL;
	std::pmr::string rep_acc, member_acc;
	AccMapping(OId rep, OId member, std::pmr::string&& member_acc, std::pmr::memory_resource& pool) :
		rep(rep),
		member(member),
		rep_acc(&pool),
		member_acc(std::move(member_acc))
	{
	}
	AccMapping(OId rep, OId member, const std::pmr::string& member_acc, std::pmr::memory_resource& pool) :
		rep(rep),
		member(member),
		rep_acc(&pool),
		member_acc(member_acc, &pool)
	{
	}
	AccMapping(std::pmr::memory_resource& pool) :
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
	friend void deserialize(File& f, AccMapping& m) {
		f.read(m.rep);
		f.read(m.member);
		f.read_c_str(m.rep_acc);
		f.read_c_str(m.member_acc);
	}
};

std::pmr::unordered_map<OId, std::pmr::string> read_mapping_table(Job& job, const Volume& vol, size_t v, std::pmr::memory_resource& pool, bool remove) {
	std::pmr::unordered_map<OId, std::pmr::string> oid2acc(&pool);
	oid2acc.reserve(vol.record_count);
	ifstream in(vol.path);
	if (!in.good())
		throw runtime_error("Error opening accessions file: " + vol.path);
	OId oid;
	std::pmr::string acc(&pool);
	while (in >> oid) {
		in >> acc;
		if (!in)
			throw runtime_error("Format error in accessions file: " + vol.path);
		if (oid2acc.emplace(oid, acc).second == false)
			throw runtime_error("Duplicate OID in accessions file: " + vol.path);
	}
	in.close();
	if (oid2acc.size() < vol.record_count)
		throw runtime_error("Accessions file does not contain all OIDs");
	if(remove)
		remove_tmp_file(vol.path);
	return oid2acc;
}

/* Rough memory estimate for one round 1 worker. It holds the accession map of a single
   volume, i.e. one hash node plus key and string header per record, and the accession
   strings themselves, for which the size of the input file is a good enough upper bound.
   On top of that come the compression buffers, one per radix, each holding an output
   buffer plus the (rather large) state of the compression stream. */
static constexpr uint64_t MAP_RECORD_MEM = 96;
static constexpr uint64_t COMPRESSION_STREAM_MEM = 320 * 1024;

static uint64_t round1_worker_mem(const VolumedFile& volumes) {
	uint64_t volume_mem = 0;
	for (const Volume& vol : volumes)
		volume_mem = std::max(volume_mem, (uint64_t)vol.record_count * MAP_RECORD_MEM + (uint64_t)file_size(vol.path.c_str()));
	return volume_mem + RADIX_COUNT * COMPRESSION_STREAM_MEM;
}

static RadixedTable output_round1(Job& job, const vector<OId>& merged, const VolumedFile& volumes) {
	const string base_dir = job.root_dir() + "output" + PATH_SEPARATOR;
	job.make_temp_dir(base_dir);
	unique_ptr<FileArray> output_files(new FileArray(base_dir, RADIX_COUNT, job.worker_id(), false));
	const int shift = std::max(bit_length(job.max_oid()) - RADIX_BITS, 0);
	const size_t volume_count = volumes.size();
	const uint64_t worker_mem = round1_worker_mem(volumes);
	const uint64_t merged_mem = (uint64_t)merged.size() * sizeof(OId);
	const uint64_t budget = job.mem_limit > merged_mem ? job.mem_limit - merged_mem : 0;
	const int64_t mem_threads = std::max<int64_t>(budget / worker_mem, 1);
	const int thread_count = (int)std::min<int64_t>(std::min<int64_t>(std::max(config.threads_, 1), (int64_t)volume_count), mem_threads);
	job.log("Building output table (round 1) threads=%i memory_estimate=%" PRIu64, thread_count, worker_mem * (uint64_t)thread_count + merged_mem);
	atomic<size_t> next(0);

	auto worker = [&](const atomic<bool>& stop) {
		std::pmr::unsynchronized_pool_resource pool;
		BufferArray buffers(*output_files, RADIX_COUNT);
		for (;;) {
			if (stop)
				return;
			const size_t v = next++;
			if (v >= volume_count)
				break;
			const Volume& vol = volumes.at(v);
			job.log("Building output table (round 1) volume %zu/%zu", v + 1, volume_count);
			const std::pmr::unordered_map<OId, std::pmr::string> oid2acc = read_mapping_table(job, vol, v, pool, false);
			remove_tmp_file(vol.path);
			for (auto it = oid2acc.cbegin(); it != oid2acc.cend(); ++it) {
				const OId oid = it->first;
				const std::pmr::string& acc = it->second;
				AccMapping m(merged[oid], oid, acc, pool);
				buffers.write(m.rep >> shift, m);
			}
			
		}
		buffers.finish();
	};

	log_rss();
	SimpleThreadPool thread_pool;
	vector<std::thread::id> workers;
	for (int i = 0; i < thread_count; ++i)
		workers.push_back(thread_pool.spawn(worker));
	thread_pool.join(workers.begin(), workers.end());
	log_rss();

	volumes.remove(true, false, true);
	return output_files->buckets(shift);
}

/* Rough memory estimate for one round 2 worker. It holds the deserialized records of a
   single bucket, i.e. a list node, the two OIDs and the two accession strings per record,
   plus the formatted output block that is queued for writing. */
static constexpr uint64_t OUTPUT_RECORD_MEM = 224;

static uint64_t round2_worker_mem(const RadixedTable& round1) {
	uint64_t records = 0;
	for (const Bucket& b : round1)
		records = std::max(records, (uint64_t)VolumedFile(b).sparse_records());
	return std::max(records * OUTPUT_RECORD_MEM, UINT64_C(1));
}

/* Takes the formatted output blocks of the round 2 workers and writes them to the output
   file in ascending bucket order, independently of the order in which the workers finish. */
struct OrderedOutput {

	OrderedOutput(std::ostream& out, Job& job, size_t bucket_count) :
		out_(out),
		job_(job),
		bucket_count_(bucket_count)
	{
	}

	/* Blocks until all preceding buckets have been written. Returns false if another worker
	   has failed in the meantime, in which case nothing is written. */
	bool write(size_t bucket, const string& block, OId clusters, size_t records) {
		std::unique_lock<std::mutex> lock(mtx_);
		cv_.wait(lock, [this, bucket]() { return aborted_ || next_ == bucket; });
		if (aborted_)
			return false;
		job_.log("Building output table (round 2) bucket %zu/%zu records=%zu", bucket + 1, bucket_count_, records);
		out_.write(block.data(), block.size());
		if (!out_.good())
			throw runtime_error("Error writing output file: " + config.output_file);
		cluster_count_ += clusters;
		++next_;
		cv_.notify_all();
		return true;
	}

	/* Releases all workers that are waiting for their turn to write. */
	void abort() {
		{
			std::lock_guard<std::mutex> lock(mtx_);
			aborted_ = true;
		}
		cv_.notify_all();
	}

	OId cluster_count() const {
		return cluster_count_;
	}

private:

	std::ostream& out_;
	Job& job_;
	const size_t bucket_count_;
	std::mutex mtx_;
	std::condition_variable cv_;
	size_t next_ = 0;
	OId cluster_count_ = 0;
	bool aborted_ = false;

};

static OId output_round2(Job& job, const vector<OId>& merged, const VolumedFile& volumes, const RadixedTable& round1, Header hdr_format) {
	ofstream out(config.output_file);
	if (!out.good())
		throw runtime_error("Error opening output file: " + config.output_file);
	if (hdr_format == Header::SIMPLE)
		out << Cluster::HEADER_LINE << endl;

	const size_t bucket_count = round1.size();
	const uint64_t worker_mem = round2_worker_mem(round1);
	const uint64_t merged_mem = (uint64_t)merged.size() * sizeof(OId);
	const uint64_t budget = job.mem_limit > merged_mem ? job.mem_limit - merged_mem : 0;
	const int64_t mem_threads = std::max<int64_t>(budget / worker_mem, 1);
	const int thread_count = (int)std::min<int64_t>(std::min<int64_t>(std::max(config.threads_, 1), (int64_t)bucket_count), mem_threads);
	job.log("Building output table (round 2) threads=%i memory_estimate=%" PRIu64, thread_count, worker_mem * (uint64_t)thread_count + merged_mem);
	log_rss();

	OrderedOutput writer(out, job, bucket_count);
	atomic<size_t> next(0);

	auto worker = [&](const atomic<bool>& stop) {
		try {
			for (;;) {
				if (stop) {
					writer.abort();
					return;
				}
				const size_t i = next++;
				if (i >= bucket_count)
					return;
				string block;
				OId clusters = 0;
				size_t records = 0;
				{
					std::pmr::unsynchronized_pool_resource pool;
					const VolumedFile f(round1[i]);
					records = (size_t)f.sparse_records();
					/* The pool is not synchronized, so the volumes of the bucket are read single
					   threaded. Parallelism comes from processing several buckets at once. */
					InputBuffer<AccMapping> data(f, pool, 1, 1);
					data.sort();
					for (auto it = data.cbegin(); it != data.cend();) {
						const auto begin = it;
						const OId rep = it->rep;
						const std::pmr::string* rep_acc = nullptr;
						while (it != data.cend() && it->rep == rep) {
							if (it->member == rep)
								rep_acc = &it->member_acc;
							++it;
						}
						if (rep_acc == nullptr)
							throw runtime_error("Missing accession mapping for representative OID " + std::to_string(rep));
						for (auto member = begin; member != it; ++member) {
							block.append(rep_acc->data(), rep_acc->size());
							block += '\t';
							block.append(member->member_acc.data(), member->member_acc.size());
							block += '\n';
						}
						++clusters;
					}
				}
				if (!writer.write(i, block, clusters, records))
					return;
			}
		}
		catch (...) {
			writer.abort();
			throw;
		}
	};

	SimpleThreadPool thread_pool;
	vector<std::thread::id> workers;
	for (int i = 0; i < thread_count; ++i)
		workers.push_back(thread_pool.spawn(worker));
	thread_pool.join(workers.begin(), workers.end());
	log_rss();

	for (const Bucket& b : round1)
		VolumedFile(b).remove();
	for (size_t v = 0; v < volumes.size(); ++v)
		remove_tmp_file(job.root_dir() + "input" + std::to_string(v) + ".tsv");
	rmdir(job.root_dir() + "output");
	return writer.cluster_count();
}

/* Writes the representative sequences to the --reps file in FASTA format. The original input
   files are read again, so the sequences are written unmasked and with their full titles. OIDs
   are assigned in input order, the same way as in make_blocks. */
static OId output_reps(Job& job, const vector<OId>& merged, const VolumedFile& input_volumes) {
	job.log("Writing representative sequences to %s", config.reps_out.c_str());
	ofstream out(config.reps_out, std::ios::out | std::ios::binary);
	if (!out.good())
		throw runtime_error("Error opening file " + config.reps_out);
	const SequenceFile::Flags flags = SequenceFile::Flags::ALL | SequenceFile::Flags::FULL_TITLES;
	const uint64_t limit = std::min<uint64_t>(job.mem_limit, SequenceFile::DEFAULT_LOAD_SIZE);
	const size_t thread_count = (size_t)std::max(config.threads_, 1);
	size_t oid = 0;
	OId count = 0;
	for (const Volume& volume : input_volumes) {
		unique_ptr<SequenceFile> file;
		try {
			file.reset(SequenceFile::auto_create({ volume.path }, flags, amino_acid_traits));
		}
		catch (FormatDetectionError& e) {
			throw runtime_error("Error opening file " + volume.path + ": " + e.what());
		}
		for (;;) {
			const unique_ptr<Block> b(file->load_seqs(limit));
			if (b->empty())
				break;
			const size_t seq_count = b->seqs().size(), oid_begin = oid;
			if (oid_begin + seq_count > merged.size())
				throw runtime_error("Input file contains more sequences than the clustering: " + volume.path);
			const size_t chunk_count = std::min(thread_count, seq_count), chunk_size = (seq_count + chunk_count - 1) / chunk_count;
			vector<TextBuffer> buffers(chunk_count);
			vector<OId> counts(chunk_count, 0);
			auto worker = [&](const atomic<bool>& stop, size_t c) {
				const size_t end = std::min(seq_count, (c + 1) * chunk_size);
				for (size_t j = c * chunk_size; j < end; ++j) {
					if (merged[oid_begin + j] != (OId)(oid_begin + j))
						continue;
					Util::Seq::format(b->seqs()[j], b->ids()[j], nullptr, buffers[c], "fasta", amino_acid_traits);
					++counts[c];
				}
			};
			SimpleThreadPool pool;
			vector<std::thread::id> workers;
			for (size_t c = 0; c < chunk_count; ++c)
				workers.push_back(pool.spawn(worker, c));
			pool.join(workers.begin(), workers.end());
			for (size_t c = 0; c < chunk_count; ++c) {
				out.write(buffers[c].data(), buffers[c].size());
				count += counts[c];
			}
			if (!out.good())
				throw runtime_error("Error writing file " + config.reps_out);
			oid += seq_count;
		}
	}
	if (oid != merged.size())
		throw runtime_error("Number of sequences in the input files does not match the clustering");
	out.close();
	if (!out)
		throw runtime_error("Error writing file " + config.reps_out);
	return count;
}

void merge(Job& job, const VolumedFile& volumes, const VolumedFile& input_volumes, Header hdr_format) {
	job.log("Merging clusterings");
	const vector<OId> merged = build_merged(job);
	OId n;
	if (config.oid_output)
		n = output_oids(job, merged);
	else {
		const RadixedTable round1 = output_round1(job, merged, volumes);
		n = output_round2(job, merged, volumes, round1, hdr_format);
	}
	job.log("Total clusters: %" PRId64, (int64_t)n);
	if (!config.reps_out.empty()) {
		const OId reps = output_reps(job, merged, input_volumes);
		if (reps != n)
			throw runtime_error("Number of representative sequences written does not match the number of clusters");
	}
}
