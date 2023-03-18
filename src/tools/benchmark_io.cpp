/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include <string>
#include <random>
#ifdef WITH_BLASTDB
#include "../data/blastdb/blastdb.h"
#endif
#include "../util/system/system.h"
#include "../util/io/output_file.h"
#include "../util/io/input_file.h"
#include "../util/log_stream.h"
#include "../data/reference.h"
#include "../util/io/input_stream_buffer.h"
#include "../basic/config.h"
#include "../util/data_structures/deque.h"
#define _REENTRANT
#include "../lib/ips4o/ips4o.hpp"
#include "../search/hit.h"

using std::vector;
using std::string;
using std::endl;

static void seed_hit_files() {
	const string file_name = "diamond_io_benchmark.tmp";
	const size_t total_count = 1000000000, query_count = 50;

	TaskTimer timer;

	if (!exists(file_name)) {
		timer.go("Writing output file");
		OutputFile out(file_name);
		std::default_random_engine generator;
		std::uniform_int_distribution<uint32_t> query(0, 2000000), seed(0, 20000), subject(1, UINT32_MAX);
		std::uniform_int_distribution<uint16_t> score(30, 1000);
		for (size_t i = 0; i < total_count / query_count; ++i) {
			out.set(Serializer::VARINT);
			out << query(generator) << seed(generator);
			out.unset(Serializer::VARINT);
			for (size_t j = 0; j < query_count; ++j) {
				out << subject(generator);
				out.write(score(generator));
			}
			out.unset(Serializer::VARINT);
			out << (uint32_t)0;
		}
		const size_t s = out.tell();
		message_stream << "Written " << (double)s / (1 << 30) << "GB. (" << s << ")" << endl;
		message_stream << "Throughput: " << (double)s / (1 << 20) / timer.seconds() << " MB/s" << endl;
		out.close();
	}

	const size_t raw_size = file_size(file_name.c_str());
	message_stream << "File size = " << raw_size << endl;
	timer.go("Reading input file");
	InputFile in(file_name, InputStreamBuffer::ASYNC);
	if (config.raw) {
		char* buf = new char[raw_size];
		in.read(buf, raw_size);
		delete[] buf;
	}
	else {
		vector<Search::Hit> out;
		out.reserve(total_count);
		auto it = std::back_inserter(out);
		size_t count = 0;
		try {
			//while (true) count += Search::Hit::read(in, it, { false });
		}
		catch (EndOfStream&) {}
		message_stream << "Read " << count << " hits." << endl;
	}
	in.close();
	timer.finish();
	message_stream << "Throughput: " << (double)raw_size / (1 << 20) / timer.seconds() << " MB/s" << endl;
}

static void load_seqs() {
	if (config.chunk_size == 0.0)
		config.chunk_size = 2.0;
	TaskTimer timer;
	timer.go("Opening the database");
	SequenceFile* db = SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NONE);
	timer.finish();
	message_stream << "Type: " << to_string(db->type()) << endl;
	Block* ref;

	while (true) {
		timer.go("Loading sequences");
		if ((ref = db->load_seqs((size_t)(config.chunk_size * 1e9)))->empty())
			return;
		size_t n = ref->seqs().letters() + ref->ids().letters();
		message_stream << "Throughput: " << (double)n / (1 << 20) / timer.milliseconds() * 1000 << " MB/s" << endl;
		timer.go("Deallocating");
		delete ref;
	}

	timer.go("Closing the database");
	db->close();
	delete db;
}

static void load_raw() {
	const size_t N = 2 * GIGABYTES;
	InputFile f(config.database);
	vector<char> buf(N);
	TaskTimer timer;
	size_t n;
	do {
		timer.go("Loading data");
		n = f.read_raw(buf.data(), N);
		timer.finish();
		message_stream << "Throughput: " << (double)n / (1 << 20) / timer.milliseconds() * 1000 << " MB/s" << endl;
	} while (n == N);
	f.close();
}

static void load_mmap() {
	static const size_t BUF = 2 * GIGABYTES;
	TaskTimer timer("Opening the database");
	SequenceFile* db = SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NONE);
	timer.finish();
	message_stream << "Type: " << to_string(db->type()) << endl;
	size_t n = db->sequence_count(), l = 0;
	vector<Letter> v, buf;
	buf.reserve(BUF);
	for (size_t i = 0; i < n; ++i) {
		db->seq_data(i, v);
		l += v.size();
		if (buf.size() + v.size() >= BUF)
			buf.clear();
		buf.insert(buf.end(), v.begin(), v.end());
		if ((i & ((1 << 20) - 1)) == 0)
			message_stream << "Throughput: " << (double)l / (1 << 20) / timer.milliseconds() * 1000 << " MB/s" << endl;
	}
	message_stream << "Throughput: " << (double)l / (1 << 20) / timer.milliseconds() * 1000 << " MB/s" << endl;
}

static void load_mmap_mt() {
	TaskTimer timer("Opening the database");
	SequenceFile* db = SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NONE);
	timer.finish();
	message_stream << "Type: " << to_string(db->type()) << endl;
	size_t n = db->sequence_count();
	std::atomic_size_t i(0);
	vector<std::thread> threads;
	for (int j = 0; j < config.threads_; ++j)
		threads.emplace_back([&i, n, db] {
		size_t k, l = 0;
		vector<Letter> v;
		while ((k = i++) < n) {
			db->seq_data(k, v);
			l += v.size();
		}
			});

	for (auto& t : threads)
		t.join();
	message_stream << "Throughput: " << (double)db->letters() / (1 << 20) / timer.milliseconds() * 1000 << " MB/s" << endl;
}

#ifdef WITH_BLASTDB
void load_blast_seqid() {
	const size_t N = 100000;
	TaskTimer timer("Opening the database");
	SequenceFile* db = SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NONE);
	timer.finish();
	message_stream << "Type: " << to_string(db->type()) << endl;
	std::mt19937 g;
	std::uniform_int_distribution<int> dist(0, (int)db->sequence_count() - 1);
	size_t n = 0;
	timer.go("Loading seqids");
	for (size_t i = 0; i < N; ++i) {
		auto l = ((BlastDB*)db)->db_->GetSeqIDs(dist(g));
		n += l.size();
		if (i % 1000 == 0)
			message_stream << i << endl;
	}
	timer.finish();
	message_stream << n << endl;
}

void load_blast_seqid_lin() {
	TaskTimer timer("Opening the database");
	SequenceFile* db = SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NONE);
	timer.finish();
	message_stream << "Type: " << to_string(db->type()) << endl;
	size_t n = 0;
	const int count = (int)db->sequence_count();
	timer.go("Loading seqids");
	for (int i = 0; i < count; ++i) {
		auto l = ((BlastDB*)db)->db_->GetSeqIDs(i);
		n += l.size();
		/*if (i % 1000 == 0)
			message_stream << i << endl;*/
	}
	timer.finish();
	message_stream << n << endl;
}

#endif

static void sort() {
	typedef uint64_t T;
	typedef Deque<T, 28> Container;
	const size_t SIZE = 1 * GIGABYTES;
	const size_t N = SIZE / sizeof(T);
	TaskTimer timer("Generating data");
	Container v;
	v.reserve(N);
	std::default_random_engine generator;
	std::uniform_int_distribution<T> r(0, std::numeric_limits<T>::max());
	for (size_t i = 0; i < N; ++i)
		v.push_back(r(generator));
	timer.go("Sorting");
	ips4o::parallel::sort(v.begin(), v.end(), std::less<T>(), config.threads_);
}

void benchmark_io() {
	if (config.type == "seedhit")
		seed_hit_files();
	else if (config.type == "loadseqs")
		load_seqs();
	else if (config.type == "loadraw")
		load_raw();
	else if (config.type == "mmap")
		load_mmap();
	else if (config.type == "mmap_mt")
		load_mmap_mt();
#ifdef WITH_BLASTDB
	else if (config.type == "blast_seqid")
		load_blast_seqid();
	else if (config.type == "blast_seqid_lin")
		load_blast_seqid_lin();
#endif
	else if (config.type == "ips4o")
		sort();
}