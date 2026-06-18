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
#include <atomic>
#include <limits>
#include <memory>
#include "multinode.h"
#include "util/data_structures/queue.h"
#include "util/parallel/simple_thread_pool.h"
#include "util/log_stream.h"
#define _REENTRANT
#include "lib/ips4o/ips4o.hpp"

using std::ifstream;
using std::unique_ptr;
using std::runtime_error;
using std::ofstream;
using std::pair;
using std::string;
using std::vector;
using std::endl;
using std::tie;
using std::atomic;

struct OutputChunk {
	static constexpr size_t POISON = std::numeric_limits<size_t>::max();
	size_t output = POISON;
	string data;
	string accs;

	bool operator==(const OutputChunk& rhs) const {
		return output == rhs.output;
	}
};

static void flush_buffer(TextBuffer& buffer, TextBuffer& acc_buffer, size_t output, vector<unique_ptr<Queue<OutputChunk>>>& queues, size_t writer_count) {
	if (buffer.size() == 0)
		return;
	OutputChunk chunk;
	chunk.output = output;
	chunk.data.assign(buffer.data(), buffer.size());
	chunk.accs.assign(acc_buffer.data(), acc_buffer.size());
	queues[output % writer_count]->enqueue(std::move(chunk));
	buffer.clear();
	acc_buffer.clear();
}

static void write_blocks(Job& job, VolumedFile& volumes, vector<unique_ptr<ofstream>>& out, vector<unique_ptr<ofstream>>& acc_out, const vector<pair<int, OId>>& block_mapping) {
	job.log("Writing length sorted blocks");
	OId oid = 0;
	uint64_t bytes_written = 0;
	const size_t writer_count = std::max<size_t>(1, std::min<size_t>(out.size(), std::max(1, config.threads_)));
	const int formatter_count = std::max(1, config.threads_);
	const size_t flush_size = 64 * 1024;
	vector<unique_ptr<Queue<OutputChunk>>> queues;
	SimpleThreadPool pool;
	vector<std::thread::id> writer_threads;
	queues.reserve(writer_count);
	for (size_t i = 0; i < writer_count; ++i)
		queues.emplace_back(new Queue<OutputChunk>(std::max<size_t>(16, formatter_count * 4), formatter_count, 1, OutputChunk()));
	auto writer = [&](const atomic<bool>& stop, size_t writer_id) {
		OutputChunk chunk;
		while (!stop.load(std::memory_order_relaxed) && queues[writer_id]->wait_and_dequeue(chunk)) {
			out[chunk.output]->write(chunk.data.data(), chunk.data.size());
			acc_out[chunk.output]->write(chunk.accs.data(), chunk.accs.size());
			bytes_written += chunk.data.size() + chunk.accs.size();
		}
	};
	for (size_t i = 0; i < writer_count; ++i)
		writer_threads.push_back(pool.spawn(writer, i));

	auto process_volume = [&](const Volume& volume, ptrdiff_t idx) {
		const SequenceFile::Flags flags = SequenceFile::Flags::ALL | SequenceFile::Flags::NEED_LETTER_COUNT;
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
			const OId oid_begin = oid;
			const size_t seq_count = b->seqs().size();
			atomic<size_t> next(0);
			vector<std::thread::id> formatter_threads;
			auto worker = [&](const atomic<bool>& stop) {
				vector<TextBuffer> buffers(out.size());
				vector<TextBuffer> acc_buffers(out.size());
				size_t j;
				while (!stop.load(std::memory_order_relaxed) && (j = next.fetch_add(1, std::memory_order_relaxed), j < seq_count)) {
					const OId seq_oid = oid_begin + (OId)j;
					const size_t block = (size_t)block_mapping[seq_oid].first;
					const string new_oid = std::to_string(block_mapping[seq_oid].second);
					Util::Seq::format(b->seqs()[j], new_oid.c_str(), nullptr, buffers[block], "fasta", amino_acid_traits);
					acc_buffers[block] << new_oid << '\t';
					acc_buffers[block] << Util::Seq::seqid(b->ids()[j]) << '\n';
					if (buffers[block].size() >= flush_size)
						flush_buffer(buffers[block], acc_buffers[block], block, queues, writer_count);
				}
				for (size_t block = 0; block < buffers.size(); ++block)
					flush_buffer(buffers[block], acc_buffers[block], block, queues, writer_count);
				};
			for (int j = 0; j < formatter_count; ++j)
				formatter_threads.push_back(pool.spawn(worker));
			try {
				pool.join(formatter_threads.begin(), formatter_threads.end());
			}
			catch (...) {
				delete b;
				for (int i = 0; i < formatter_count; ++i)
					for (auto& q : queues)
						q->close();
				pool.join(writer_threads.begin(), writer_threads.end());
				throw;
			}
			oid += (OId)seq_count;
			delete b;
		}
		};

	TaskTimer timer;
	for (VolumedFile::const_iterator i = volumes.begin(); i != volumes.end(); ++i) {
		process_volume(*i, i - volumes.begin());
	}
	const int64_t t = timer.microseconds();
	job.log("Wrote %zu bytes to disk at %.2f MB/s", bytes_written, (double)bytes_written / MEGABYTES / (t / 1e6));
	for (int i = 0; i < formatter_count; ++i)
		for (auto& q : queues)
			q->close();
	pool.join(writer_threads.begin(), writer_threads.end());
	for (const auto& file : out)
		if (!*file)
			throw runtime_error("Error writing length sorted block");
	for (const auto& file : acc_out)
		if (!*file)
			throw runtime_error("Error writing length sorted block");
}

string len_sort(Job& job, VolumedFile& volumes) {
	Atomic lock(job.root_dir() + "lensort_lock", job), done(job.root_dir() + "lensort_done", job);
	const bool first_round_linear = job.is_linear_round();
	const string input_parts = job.root_dir() + "input.tsv", input_letters = job.root_dir() + "input_letters.txt";
	if (lock.fetch_add() == 0) {
		job.log("Memory limit = %" PRIu64, job.mem_limit);
		vector<pair<Loc, OId>> lengths;
		uint64_t letters = 0;
		for (vector<Volume>::const_iterator i = volumes.begin(); i != volumes.end(); ++i) {
			job.log("Reading volume %td/%zu", i - volumes.begin() + 1, volumes.size());
			unique_ptr<SequenceFile> volume_file;
			try {
				volume_file.reset(SequenceFile::auto_create({ i->path }, SequenceFile::Flags::NEED_LETTER_COUNT | SequenceFile::Flags::NEED_LENGTH_LOOKUP, amino_acid_traits));
			} catch(FormatDetectionError& e) {
				throw runtime_error("Error opening file " + i->path + ": " + e.what());
			}
			string msg = volume_file->open_stats();
			if (volume_file->type() == SequenceFile::Type::DMND)
				job.log("Warning: using legacy loader on .dmnd files");
			else {
				msg.pop_back();
				job.log(msg.c_str());
			}
			letters += volume_file->letters().value();
			const OId n = volume_file->sequence_count().value();
			for (OId i = 0; i < n; ++i)
				lengths.emplace_back(volume_file->seq_length(i), lengths.size());
		}
		ofstream letters_out(input_letters);
		letters_out << letters << endl;
		letters_out << lengths.size() << endl;
		if (!letters_out)
			throw runtime_error("Error writing file " + input_letters);
		letters_out.close();
		job.log("Computing blocks");
		std::ostringstream ss;
		//volumes.set_max_oid(lengths.size() - 1);
		ss << "Sequences in database = " << lengths.size() << endl;
		ss << "Letters in database = " << letters << endl;
		ss << "Database blocks:" << endl;
		
		double block_gb;
		int index_chunks;
		if (!first_round_linear) {
			block_gb = 1e6;
			index_chunks = 1;
		} else
			tie(block_gb, index_chunks) = ::block_size(job.mem_limit, letters, Sensitivity::FAST, true, config.threads_); // TODO take cluster steps into account here
		const uint64_t block_size = gb_to_bytes(block_gb);
		job.log("Block size = %" PRIu64 ", index chunks = %d", block_size, index_chunks);
		ips4o::parallel::sort(lengths.begin(), lengths.end(), std::greater<pair<Loc, OId>>(), config.threads_);
		ofstream idx(input_parts);
		vector<unique_ptr<ofstream>> out;
		vector<unique_ptr<ofstream>> acc_out;
		vector<pair<int, OId>> block_mapping(lengths.size());
		int block = 0;
		OId new_oid = 0;
		for (vector<pair<Loc, OId>>::const_iterator i = lengths.begin(); i != lengths.end();) {
			uint64_t block_letters = 0, seqs = 0;
			while (i != lengths.end() && block_letters + i->first <= block_size) {
				block_letters += i->first;
				block_mapping[i->second] = { block, new_oid++ };
				++i;
				++seqs;
			}
			ss << seqs << '\t' << block_letters << endl;
			const string block_idx = std::to_string(block);
			const string name = job.root_dir() + "input" + block_idx + ".faa";
			out.emplace_back(new ofstream(name));
			acc_out.emplace_back(new ofstream(job.root_dir() + "input" + block_idx + ".tsv"));
			idx << name << '\t' << seqs << endl;
			++block;
		}
		job.log(ss.str().c_str());
		write_blocks(job, volumes, out, acc_out, block_mapping);
		done.fetch_add();
	}
	else
		done.await(1);
	uint64_t letters, seq_count;
	ifstream input_letters_file(input_letters);
	input_letters_file >> letters >> seq_count;
	if(!input_letters_file)
		throw runtime_error("Error opening file " + input_letters);
	volumes.set_letter_count(letters);
	job.register_sync_file(input_letters);
	if (first_round_linear) {
		double block_gb;
		int index_chunks;
		tie(block_gb, index_chunks) = ::block_size(job.mem_limit, letters, Sensitivity::FAST, true, config.threads_); // TODO take cluster steps into account here
		config.lowmem_ = index_chunks;
	}
	job.set_max_oid(seq_count - 1);
	return input_parts;
}