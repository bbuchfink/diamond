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
#include <algorithm>
#include <atomic>
#include <exception>
#include <limits>
#include <memory>
#include <unordered_map>
#include "masking/masking.h"
#include "multinode.h"
#include "util/data_structures/queue.h"
#include "util/parallel/simple_thread_pool.h"
#include "util/log_stream.h"

using std::ifstream;
using std::unique_ptr;
using std::runtime_error;
using std::ofstream;
using std::string;
using std::vector;
using std::pair;
using std::endl;
using std::tie;
using std::atomic;

struct OutputPart {
	size_t file;
	string data;
	string accs;
};

struct OutputChunk {
	static constexpr size_t POISON = std::numeric_limits<size_t>::max();
	size_t order = POISON;
	vector<OutputPart> parts;

	bool operator==(const OutputChunk& rhs) const {
		return order == rhs.order;
	}
};

static void write_blocks(Job& job, VolumedFile& volumes, vector<unique_ptr<ofstream>>& out, vector<unique_ptr<ofstream>>& acc_out, const vector<int>& block_mapping) {
	job.log("Writing length sorted blocks");
	/* The input is masked here, once, and the masked letters are written to the minichunks
	   as X. All later stages of the workflow (seed counting, seed index building, alignment)
	   read the minichunks and therefore need no masking of their own. */
	const MaskingAlgo masking = input_masking_algo();
	job.log("Input masking: %s", to_string(masking).c_str());
	OId oid = 0;
	atomic<uint64_t> bytes_written{ 0 };
	const size_t writer_count = std::max<size_t>(1, std::min<size_t>(out.size(), std::max(1, config.threads_)));
	const int formatter_count = std::max(1, config.threads_);
	SimpleThreadPool pool;

	auto process_volume = [&](const Volume& volume, ptrdiff_t idx) {
		const SequenceFile::Flags flags = SequenceFile::Flags::ALL;
		unique_ptr<SequenceFile> file;
		try {
			file.reset(SequenceFile::auto_create({ volume.path }, flags, amino_acid_traits));
		} catch(FormatDetectionError& e) {
			throw runtime_error("Error opening file " + volume.path + ": " + e.what());
		}
		const uint64_t limit = std::min<uint64_t>(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)), SequenceFile::DEFAULT_LOAD_SIZE);
		uint64_t bytes = 0, ms = 0;
		for (;;) {
			TaskTimer timer;
			Block* b = file->load_seqs(limit);
			ms += timer.microseconds();
			bytes += b->raw_bytes();
			if (b->empty()) {
				delete b;
				break;
			}
			if (masking != MaskingAlgo::NONE) {
				timer.go("Masking input block");
				const MaskingStat stats = mask_seqs(b->seqs(), Masking::get(), true, masking);
				timer.finish();
				stats.print(*message_stream);
			}
			const OId oid_begin = oid;
			const size_t seq_count = b->seqs().size();			
			const size_t target_chunks = std::max<size_t>(1, (size_t)formatter_count * 16);
			const size_t chunk_size = std::max<size_t>(1, std::min<size_t>(seq_count / target_chunks + 1, (size_t)8192));
			const size_t chunk_count = (seq_count + chunk_size - 1) / chunk_size;
			
			vector<unique_ptr<Queue<OutputChunk>>> queues;
			queues.reserve(writer_count);
			for (size_t i = 0; i < writer_count; ++i)
				queues.emplace_back(new Queue<OutputChunk>(std::max<size_t>(16, formatter_count * 4), formatter_count, 1, OutputChunk()));

			auto writer = [&](const atomic<bool>& stop, size_t writer_id) {
				std::unordered_map<size_t, OutputChunk> pending;
				size_t next_expected = 0;
				OutputChunk chunk;
				while (!stop.load(std::memory_order_relaxed) && queues[writer_id]->wait_and_dequeue(chunk)) {
					pending.emplace(chunk.order, std::move(chunk));
					auto it = pending.find(next_expected);
					while (it != pending.end()) {
						for (OutputPart& p : it->second.parts) {
							out[p.file]->write(p.data.data(), p.data.size());
							acc_out[p.file]->write(p.accs.data(), p.accs.size());
							bytes_written.fetch_add(p.data.size() + p.accs.size(), std::memory_order_relaxed);
						}
						pending.erase(it);
						++next_expected;
						it = pending.find(next_expected);
					}
				}
			};

			atomic<size_t> next(0);
			auto worker = [&](const atomic<bool>& stop) {
				vector<TextBuffer> buffers(out.size());
				vector<TextBuffer> acc_buffers(out.size());
				size_t c;
				while (!stop.load(std::memory_order_relaxed) && (c = next.fetch_add(1, std::memory_order_relaxed), c < chunk_count)) {
					const size_t j_begin = c * chunk_size;
					const size_t j_end = std::min(j_begin + chunk_size, seq_count);
					for (size_t j = j_begin; j < j_end; ++j) {
						const OId seq_oid = oid_begin + (OId)j;
						const size_t block = (size_t)block_mapping[seq_oid];
						const string oid_str = std::to_string(seq_oid);
						Util::Seq::format(b->seqs()[j], oid_str.c_str(), nullptr, buffers[block], "fasta", amino_acid_traits);
						acc_buffers[block] << oid_str << '\t';
						acc_buffers[block] << Util::Seq::seqid(b->ids()[j]) << '\n';
					}					
					vector<OutputChunk> bundles(writer_count);
					for (size_t w = 0; w < writer_count; ++w)
						bundles[w].order = c;
					for (size_t block = 0; block < buffers.size(); ++block) {
						if (buffers[block].size() == 0)
							continue;
						OutputPart part;
						part.file = block;
						part.data.assign(buffers[block].data(), buffers[block].size());
						part.accs.assign(acc_buffers[block].data(), acc_buffers[block].size());
						bundles[block % writer_count].parts.push_back(std::move(part));
						buffers[block].clear();
						acc_buffers[block].clear();
					}
					for (size_t w = 0; w < writer_count; ++w)
						queues[w]->enqueue(std::move(bundles[w]));
				}
				};

			vector<std::thread::id> writer_threads;
			for (size_t i = 0; i < writer_count; ++i)
				writer_threads.push_back(pool.spawn(writer, i));
			vector<std::thread::id> formatter_threads;
			for (int j = 0; j < formatter_count; ++j)
				formatter_threads.push_back(pool.spawn(worker));
			try {
				pool.join(formatter_threads.begin(), formatter_threads.end());
			}
			catch (...) {
				for (int i = 0; i < formatter_count; ++i)
					for (auto& q : queues)
						q->close();
				pool.join(writer_threads.begin(), writer_threads.end());
				delete b;
				throw;
			}			
			for (int i = 0; i < formatter_count; ++i)
				for (auto& q : queues)
					q->close();
			pool.join(writer_threads.begin(), writer_threads.end());
			oid += (OId)seq_count;
			delete b;
			*message_stream << "Processed " << oid << " / " << block_mapping.size() << " sequences (" << std::setprecision(2) << std::fixed
				<< (double)oid / std::max<size_t>(block_mapping.size(), 1) * 100 << "%)" << endl;
		}
		};

	TaskTimer timer;
	for (VolumedFile::const_iterator i = volumes.begin(); i != volumes.end(); ++i) {
		process_volume(*i, i - volumes.begin());
	}
	const int64_t t = timer.microseconds();
	const uint64_t total_bytes = bytes_written.load();
	job.log("Wrote %zu bytes to disk at %.2f MB/s", total_bytes, (double)total_bytes / MEGABYTES / (t / 1e6));
	for (const auto& file : out)
		if (!*file)
			throw runtime_error("Error writing length sorted block");
	for (const auto& file : acc_out)
		if (!*file)
			throw runtime_error("Error writing length sorted block");
}

pair<string, string> len_sort(Job& job, VolumedFile& volumes) {
	Atomic lock(job.root_dir() + "lensort_lock", job), done(job.root_dir() + "lensort_done", job);
	const bool first_round_linear = job.is_linear_round();
	const string input_letters = job.root_dir() + "input_letters.txt", out_dir = job.root_dir() + "input_minichunks" + PATH_SEPARATOR,
		out_seqs = out_dir + "seqs.tsv", out_accs = out_dir + "accs.tsv";
	if (lock.fetch_add() == 0) {
		vector<unique_ptr<ofstream>> out;
		vector<unique_ptr<ofstream>> acc_out;
		vector<int> block_mapping = make_blocks(job, volumes, out, acc_out, out_seqs, out_accs);
		write_blocks(job, volumes, out, acc_out, block_mapping);
		done.fetch_add();
		job.finish_step();
	}
	else
		done.await(1);
	uint64_t letters, seq_count, minichunk_count;
	ifstream input_letters_file(input_letters);
	input_letters_file >> letters >> seq_count >> minichunk_count;
	if(!input_letters_file)
		throw runtime_error("Error opening file " + input_letters);
	volumes.set_letter_count(letters);
	job.register_sync_file(input_letters);
	if (first_round_linear && config.mutual_cover.present()) {
		double block_gb;
		int index_chunks;
		tie(block_gb, index_chunks) = ::block_size(job.mem_limit, letters, Sensitivity::FAST, true, config.threads_, config.mutual_cover.present()); // TODO take cluster steps into account here
		config.lowmem_ = index_chunks;
	}
	job.set_max_oid(seq_count - 1);
	return std::make_pair(out_seqs, out_accs);
}