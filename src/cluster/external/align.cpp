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

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
#include <unordered_map>
#include "basic/config.h"
#include "external.h"
#include "util/parallel/atomic.h"
#include "util/util.h"
#include "dp/dp.h"
#include "util/parallel/thread_pool.h"
#include "util/hash_function.h"

using std::mutex;
using std::thread;
using std::vector;
using std::string;
using std::unique_ptr;
using std::atomic;
using std::unordered_map;
using std::endl;
using std::lock_guard;
using std::list;
using std::pair;
using std::array;

namespace Cluster {

static void align_rep(ThreadPool& tp, const ChunkSeqs& chunk_seqs, vector<PairEntryShort>::const_iterator begin, vector<PairEntryShort>::const_iterator end, int64_t db_size, BufferArray& out) {
	const int64_t rep_oid = begin->rep_oid, shift = bit_length(db_size - 1) - RADIX_BITS;
	const Sequence rep = chunk_seqs[rep_oid];
	DP::Targets targets;
	Statistics stats;
	Loc max_len = 0;
	for (auto i = begin; i < end; ++i) {
		const Sequence member = chunk_seqs[i->member_oid];
		const int bin = DP::BandedSwipe::bin(HspValues::COORDS, rep.length(), 0, 0, (int64_t)rep.length() * (int64_t)member.length(), 0, 0);
		targets[bin].emplace_back(member, member.length(), BlockId(i - begin));
		max_len = std::max(max_len, member.length());
	}
	DP::Params params{ rep, nullptr, Frame(0), rep.length(), nullptr, DP::Flags::FULL_MATRIX, false, max_len,
		-1, HspValues::COORDS, stats, &tp };
	const list<Hsp> hsps = DP::BandedSwipe::swipe(targets, params);
	const bool unid = !config.mutual_cover.present();
	for (const Hsp& h : hsps) {
		const int64_t member_oid = begin[h.swipe_target].member_oid;
		const Sequence member = chunk_seqs[member_oid];
		if (h.approx_id_percent(rep, member) < config.approx_min_id.get(0))
			continue;
		if (unid) {
			if (h.subject_cover_percent(member.length()) >= config.member_cover.get(80))
				out.write(MurmurHash()(member_oid) & (RADIX_COUNT - 1), Edge(rep_oid, member_oid, rep.length(), member.length()));
			if (h.query_cover_percent(rep.length()) >= config.member_cover.get(80))
				out.write(MurmurHash()(rep_oid) & (RADIX_COUNT - 1), Edge(member_oid, rep_oid, member.length(), rep.length()));
		}
		else if (h.subject_cover_percent(member.length()) >= config.mutual_cover.get_present() && h.query_cover_percent(rep.length()) >= config.mutual_cover.get_present()) {
			int64_t oid1 = rep_oid;
			int64_t oid2 = member_oid;
			Loc len1 = rep.length();
			Loc len2 = member.length();
			if (oid1 > oid2) {
				std::swap(oid1, oid2);
				std::swap(len1, len2);
			}
			out.write(int(MurmurHash()(oid1) & (RADIX_COUNT - 1)), Edge(oid1, oid2, len1, len2));
		}
	}
}

vector<string> align(Job& job, int chunk_count, int64_t db_size) {
	const std::string chunks_path = job.base_dir() + PATH_SEPARATOR + "chunks" + PATH_SEPARATOR,
		queue_path = chunks_path + "align_queue",
		aln_path = job.base_dir() + PATH_SEPARATOR + "alignments";
	score_matrix.set_db_letters(1000000000);
	mkdir(aln_path);
	unique_ptr<FileArray> output_files(new FileArray(aln_path, RADIX_COUNT, job.worker_id()));
	Atomic queue(queue_path);
	int chunk, chunks_processed = 0;
	while (true) {
		chunk = (int)queue.fetch_add();
		if (chunk >= chunk_count)
			break;
		TaskTimer timer("Reading sequence files");
		const string chunk_path = chunks_path + std::to_string(chunk) + PATH_SEPARATOR;
		unique_ptr<ChunkSeqs> chunk_seqs(new ChunkSeqs(chunk_path));
		timer.finish();
		job.log("Computing alignments. Chunk=%lli/%lli Volumes=%lli Sequences=%s Letters=%s", chunk + 1, chunk_count, chunk_seqs->volumes(), Util::String::format(chunk_seqs->oids()).c_str(), Util::String::format(chunk_seqs->letters()).c_str());
		timer.go("Computing alignments");
		InputFile pairs_file(chunk_path + "pairs", InputFile::NO_AUTODETECT);
		int64_t pairs_processed = 0;
		ThreadPool tp;
		ThreadPool::TaskSet task_set(tp, 1);
		std::function<void()> align_worker = [&]() {
			thread_local vector<PairEntryShort> buf;
			thread_local BufferArray edge_buf(*output_files, RADIX_COUNT);
			size_t size;
			if (pairs_file.read(&size, 1) == 0)
				return;
			buf.clear();
			buf.resize(size);
			pairs_file.read(buf.data(), size);
			pairs_processed += size;
			task_set.enqueue(align_worker);
			auto it = merge_keys(buf.begin(), buf.end(), PairEntryShort::Key());
			while (it.good()) {
				align_rep(tp, *chunk_seqs, it.begin(), it.end(), db_size, edge_buf);
				++it;
			}
			};
		task_set.enqueue(align_worker);
		tp.run(config.threads_, true, &task_set);
		tp.join();
		timer.finish();
		log_stream << "pairs=" << pairs_processed << endl;
		timer.go("Deallocating memory");
		pairs_file.close();
		remove(pairs_file.file_name.c_str());
		chunk_seqs.reset();
		++chunks_processed;
		timer.finish();
	}
	TaskTimer timer("Closing the output files");
	const vector<string> out = output_files->buckets();
	output_files.reset();
	timer.go("Waiting for other workers");
	Atomic finished(aln_path + PATH_SEPARATOR + "finished");
	finished.fetch_add(chunks_processed);
	finished.await(chunk_count);
	return out;
	//const auto ccs = cc.largest_cc();
	//log_stream << "Connected components = " << ccs.second << endl;
	//log_stream << "Largest connected component = " << ccs.first << endl;
	//TaskTimer timer("Allocating buffers");
	//vector<Edge<int64_t>> edges(edge_count);
	//vector<Edge> edges(edge_count);
	/*timer.go("Loading edges");
	output_file->close();
	InputFile f(aln_path + "worker_0");
	f.read(edges.data(), edge_count);
	f.close();
	timer.go("Sorting edges");*/
	//FlatArray<Edge<int64_t>> edge_array = make_flat_array_dense(std::move(edges), db_size, config.threads_, Edge<int64_t>::GetKey());
	//ips4o::parallel::sort(edges.begin(), edges.end(), std::less<Edge>(), config.threads_);
	//timer.finish();
	//vector<int64_t> oid2rep = Util::Algo::greedy_vertex_cover<int64_t>(edge_array, nullptr, true, false, 0);
	/*vector<int64_t> rep(db_size);
	std::iota(rep.begin(), rep.end(), 0);
	auto it = merge_keys(edges.begin(), edges.end(), Edge::Member());
	while (it.good()) {
		if (it.begin()->member_len < it.begin()->rep_len || (it.begin()->member_len == it.begin()->rep_len && it.begin()->member_oid > it.begin()->rep_oid))
			rep[it.begin()->member_oid] = it.begin()->rep_oid;
		++it;
	}
	int64_t changed = 0;
	do {
		changed = 0;
		for (int64_t i = 0; i < db_size; ++i)
			if (rep[rep[i]] != rep[i]) {
				rep[i] = rep[rep[i]];
				++changed;
			}
	} while (changed > 0);
	int64_t n = 0;
	for (int64_t i = 0; i < db_size; ++i)
		if (rep[i] == i)
			++n;
	log_stream << n << endl;*/
}

}