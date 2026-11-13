/****
DIAMOND protein aligner
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

#include <fstream>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <memory>
#include <queue>
#define _REENTRANT
#include "util/algo/degree_partition.h"
#include "lib/ips4o/ips4o.hpp"
#include "basic/config.h"
#include "util/log_stream.h"
#include "util/string/tokenizer.h"
#include "util/string/string.h"
#include "util/algo/algo.h"
#include "util/system/system.h"
#include "cluster/cluster.h"
#include "cluster/file_array.h"
#include "cluster/input_buffer.h"
#include "util/memory/memory_resource.h"
#include "util/parallel/simple_thread_pool.h"
#include "tools.h"
#include "data/fasta/parser.h"

struct Edge {
	static constexpr bool POD = true;
	Edge()
	{
	}
	using Key = OId;
	Edge(Key node1, Key node2, double weight) :
		node1(node1),
		node2(node2),
		weight(weight)
	{
	}
	bool operator<(const Edge& e) const noexcept {
		if (node1 != e.node1)
			return node1 < e.node1;
		if (node2 != e.node2)
			return node2 < e.node2;
		return weight > e.weight;
	}
	struct GetKey {
		OId operator()(const Edge& e) const {
			return e.node1;
		}
	};
	OId node1, node2;
	double weight;
};

inline void serialize(const Edge& e, CompressedBuffer& buf) {
	buf.write(e.node1);
	buf.write(e.node2);
	buf.write(e.weight);
}

using std::unique_ptr;
using std::ofstream;
using std::mutex;
using std::lock_guard;
using std::endl;
using std::string;
using std::atomic;
using std::vector;
using std::unordered_map;
using std::numeric_limits;
using std::runtime_error;
using std::function;
using std::priority_queue;
using std::numeric_limits;
using std::deque;
using std::pair;
using Acc = string;

namespace GVC {

static const bool merge_recursive = true;
static const uint64_t RADIX_COUNT = 256;
// Resolution of the edge count histogram of edge_pass_zero: at most 2^BIN_BITS key bins.
static const int BIN_BITS = 22;
// Share of the memory limit that one bucket of edges loaded by edge_pass_two may use.
static const double BUCKET_MEMORY_FRACTION = 0.5;

struct PotentialRep {
	OId oid, degree;
	vector<pair<OId, double>> members;
	PotentialRep(OId rep):
		oid(rep),
		degree(0)
	{
	}
	bool operator<(const PotentialRep& r) const {
		return degree < r.degree || (degree == r.degree && oid < r.oid);
	}
	void set_degree(const vector<OId>& clustering) {
		degree = 0;
		for (const auto& m : members) {
			if (clustering[m.first] == numeric_limits<OId>::max())
				++degree;
		}
	}
};

using RepQueue = priority_queue<PotentialRep, deque<PotentialRep>>;

static void greedy_vertex_cover(vector<OId>& clustering, vector<double>& weights, RepQueue& reps, uint64_t& edges_queued, OId next_degree) {
	while (!reps.empty()) {
		PotentialRep r = reps.top();
		reps.pop();
		edges_queued -= r.members.size();
		if (clustering[r.oid] != numeric_limits<OId>::max())
			continue;
		r.set_degree(clustering);
		if (!reps.empty() && r.degree < reps.top().degree) {
			reps.push(r);
			edges_queued += r.members.size();
			continue;
		}
		if (r.degree < next_degree) {
			reps.push(r);
			edges_queued += r.members.size();
			return;
		}
		clustering[r.oid] = r.oid;
		for (const auto& m : r.members) {
			if (clustering[m.first] == numeric_limits<OId>::max()
				|| (weights[m.first] < m.second && !config.no_gvc_reassign && clustering[m.first] != m.first)
				|| (merge_recursive && clustering[m.first] == m.first)) {
				clustering[m.first] = r.oid;
				weights[m.first] = m.second;
			}
		}
	}
}

struct EdgeInput {
	OId max_oid;
	bool triplets, symmetric;
	double cov;
	const unordered_map<Acc, OId>* acc2oid;
	const VolumedFile* edges;
};

// Parses the input and calls emit(node1, node2, weight) for every directed edge, concurrently from multiple threads.
// Returns the number of input lines.
template<typename Emit>
static uint64_t read_edges(const EdgeInput& in, Emit emit) {
	atomic<uint64_t> total_lines(0);
	function<void(int64_t, const char*, const char*)> fn2([&](int64_t, const char* begin, const char* end) {
		size_t line_count = 0;
		Util::String::LineIterator it(begin, end);
		string query, target;
		OId q, t;
		double qcov, tcov, weight;
		while (it.good()) {
			++line_count;
			string line = *it;
			if (ends_with(line, "\r"))
				line.pop_back();
			Util::String::Tokenizer<Util::String::CharDelimiter> tok(line, Util::String::CharDelimiter('\t'));
			if (in.acc2oid)
				tok >> query >> target;
			else {
				tok >> q >> t;
				if (q > in.max_oid || t > in.max_oid)
					throw runtime_error("Numeric vertex ID exceeds --max-id");
			}
			if (!in.triplets)
				tok >> qcov >> tcov;
			tok >> weight;
			if (in.triplets || tcov >= in.cov || qcov >= in.cov) {
				if (in.acc2oid) {
					q = in.acc2oid->at(query);
					t = in.acc2oid->at(target);
				}
				if (q == t) {
					++it;
					continue;
				}
				if (in.triplets) {
					emit(t, q, weight);
					if (in.symmetric)
						emit(q, t, weight);
				}
				else {
					if (tcov >= in.cov)
						emit(q, t, weight);
					if (qcov >= in.cov)
						emit(t, q, weight);
				}
			}
			++it;
		}
		total_lines.fetch_add(line_count, std::memory_order_relaxed);
		});
	if (in.edges) {
		for (const Volume& volume : *in.edges) {
			File f(volume.path, "rb");
			f.read_text_mt(INT64_MAX, config.threads_, fn2);
		}
	}
	else {
		File f(config.edges, "rb");
		f.read_text_mt(INT64_MAX, config.threads_, fn2);
	}
	return total_lines.load(std::memory_order_relaxed);
}

// Contiguous ranges of node1 keys, sized so that the edges of one range fit into the memory budget of edge_pass_two.
struct EdgePartition {
	size_t size() const {
		return key_bounds.size() - 1;
	}
	uint32_t bucket_of(OId key) const {
		return bucket[key >> shift];
	}
	int shift;
	vector<uint32_t> bucket;     // key bin -> bucket
	vector<uint64_t> key_bounds; // bucket i holds the keys [key_bounds[i], key_bounds[i + 1])
};

// Counts the edges per key bin and cuts the key range greedily into the fewest buckets that fit into the memory budget.
static EdgePartition edge_pass_zero(const EdgeInput& in) {
	TaskTimer timer("Counting input edges");
	EdgePartition p;
	p.shift = std::max(bit_length(in.max_oid) - BIN_BITS, 0);
	const size_t bins = size_t(in.max_oid >> p.shift) + 1;
	vector<atomic<uint64_t>> counts(bins);
	for (auto& c : counts)
		c.store(0, std::memory_order_relaxed);
	read_edges(in, [&](OId n1, OId, double) {
		counts[n1 >> p.shift].fetch_add(1, std::memory_order_relaxed);
		});
	timer.finish();

	const uint64_t budget = uint64_t(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)) * BUCKET_MEMORY_FRACTION);
	const uint64_t max_records = std::max<uint64_t>(budget / sizeof(Edge), 1);
	p.bucket.resize(bins);
	p.key_bounds.push_back(0);
	uint64_t records = 0, largest = 0, total = 0;
	for (size_t i = 0; i < bins; ++i) {
		const uint64_t c = counts[i].load(std::memory_order_relaxed);
		if (records > 0 && records + c > max_records) {
			p.key_bounds.push_back(uint64_t(i) << p.shift);
			largest = std::max(largest, records);
			records = 0;
		}
		records += c;
		total += c;
		p.bucket[i] = uint32_t(p.key_bounds.size() - 1);
	}
	largest = std::max(largest, records);
	p.key_bounds.push_back(uint64_t(in.max_oid) + 1);
	*message_stream << "#Counted edges: " << total << endl;
	*message_stream << "Edge buckets: " << p.size() << " (largest=" << largest << " edges, limit=" << max_records << " edges)" << endl;
	if (largest > max_records)
		*message_stream << "Warning: the edges of one key range exceed the memory budget of " << budget << " bytes." << endl;
	return p;
}

static RadixedTable edge_pass_one(const string& base_dir, const EdgeInput& in, const EdgePartition& p) {
	mkdir(base_dir);
	FileArray file_array(base_dir, (int)p.size(), 0, true);
	// A buffer is flushed at BufferArray::BUF_SIZE and grows in steps of CompressedBuffer::BUF_SIZE.
	const uint64_t buffer_mem = CompressedBuffer::stream_memory() + BufferArray::BUF_SIZE + CompressedBuffer::BUF_SIZE;
	*message_stream << "Edge write buffer memory: " << (uint64_t)config.threads_ * p.size() * buffer_mem << " bytes" << endl;

	TaskTimer timer("Reading input lines");
	const uint64_t total_lines = read_edges(in, [&](OId n1, OId n2, double w) {
		thread_local BufferArray buffers(file_array, (int)p.size());
		Edge edge(n1, n2, w);
		buffers.write(p.bucket_of(n1), &edge, 1, 1);
		});
	timer.finish();
	log_rss();
	*message_stream << "#Input lines: " << total_lines << endl;
	*message_stream << "#Input edges: " << file_array.records_total() << endl;
	file_array.close();
	return file_array.buckets(p.key_bounds);
}

static OId node_degree(vector<Edge>::const_iterator begin, vector<Edge>::const_iterator end) {
	OId count = 0;
	for (auto j = begin; j != end; ++j) {
		if (j->node2 != j->node1 && (j == begin || j->node2 != (j - 1)->node2))
			++count;
	}
	return count;
}

static DegreePartition edge_pass_two(const RadixedTable& rep_sorted) {
	TaskTimer timer;
	unordered_map<OId, OId> degrees;
	for (RadixedTable::const_iterator it = rep_sorted.begin(); it != rep_sorted.end(); ++it) {
		VolumedFile f(*it);
		InputBuffer<Edge> data(f);
		*message_stream << "Finding neighbor counts bucket " << it - rep_sorted.begin() + 1 << "/" << rep_sorted.size() << " records=" << data.size() << endl;
#ifdef NDEBUG
		ips4o::parallel::sort(data.begin(), data.end());
#else
		std::sort(data.begin(), data.end());
#endif
		const auto parts = Util::Algo::partition_table(data.cbegin(), data.cend(), config.threads_, Edge::GetKey());
		if (parts.empty())
			continue;
		vector<unordered_map<OId, OId>> local_degrees(parts.size() - 1);
		SimpleThreadPool pool;
		for (size_t part = 0; part < local_degrees.size(); ++part) {
			pool.spawn([&, part](const atomic<bool>& stop) {
				auto i = merge_keys(parts[part], parts[part + 1], Edge::GetKey());
				while (i.good() && !stop.load(std::memory_order_relaxed)) {
					const OId d = node_degree(i.begin(), i.end());
					local_degrees[part][d] += d;
					++i;
				}
				});
		}
		pool.join_all();
		for (const auto& local : local_degrees)
			for (const auto& degree : local)
				degrees[degree.first] += degree.second;
	}
	// With a single bucket, edge_pass_four holds all edges and the members of all potential reps at the same time,
	// so its memory use is known in advance. Otherwise the amount of edges buffered in the queue is unknown.
	uint64_t total_edges = 0;
	for (const auto& d : degrees)
		total_edges += d.second;
	const uint64_t memory_limit = Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT));
	const bool single_bucket = total_edges <= memory_limit / (sizeof(Edge) + sizeof(pair<OId, double>));
	DegreePartition r(degrees, single_bucket ? 1 : RADIX_COUNT);
	*message_stream << "Degree buckets: " << r.size() << " (edges=" << total_edges << ')' << endl;
	*message_stream << "Time (finding neighbor counts): " << timer.seconds() << 's' << endl;
	return r;
}

static RadixedTable edge_pass_three(const RadixedTable& rep_sorted, const DegreePartition& p, Cfg& cfg) {
	TaskTimer timer;
	const string base_dir = cfg.tmp_dir + "degree_sorted";
	mkdir(base_dir);
	FileArray file_array(base_dir, p.size(), 0, true);
	for (RadixedTable::const_iterator it = rep_sorted.begin(); it != rep_sorted.end(); ++it) {
		VolumedFile f(*it);
		InputBuffer<Edge> data(f);
		*message_stream << "Writing degree sorted edges " << it - rep_sorted.begin() + 1 << "/" << rep_sorted.size() << " records=" << data.size() << endl;
#ifdef NDEBUG
		ips4o::parallel::sort(data.begin(), data.end());
#else
		std::sort(data.begin(), data.end());
#endif
		const auto parts = Util::Algo::partition_table(data.cbegin(), data.cend(), config.threads_, Edge::GetKey());
		SimpleThreadPool pool;
		for (size_t part = 0; part + 1 < parts.size(); ++part) {
			pool.spawn([&, part](const atomic<bool>& stop) {
				BufferArray buffers(file_array, p.size());
				auto i = merge_keys(parts[part], parts[part + 1], Edge::GetKey());
				while (i.good() && !stop.load(std::memory_order_relaxed)) {
					const int bucket = p.bucket_index(node_degree(i.begin(), i.end()));
					for (auto j = i.begin(); j != i.end(); ++j) {
						if (j->node2 != j->node1 && (j == i.begin() || j->node2 != (j - 1)->node2))
							buffers.write(bucket, *j);
					}
					++i;
				}
				buffers.finish();
				});
		}
		pool.join_all();
		f.remove();
	}
	file_array.close();
	rmdir(cfg.tmp_dir + "rep_sorted");
	*message_stream << "Time (writing degree sorted edges): " << timer.seconds() << 's' << endl;
	return file_array.buckets(p);
}

static void edge_pass_four(const RadixedTable& degree_sorted, OId db_size, Cfg& cfg) {
	TaskTimer timer;
	cfg.clustering = std::make_unique<vector<OId>>(db_size, numeric_limits<OId>::max());
	vector<double> weights(db_size);
	RepQueue queue;
	uint64_t edges_queued = 0;
	for (int i = degree_sorted.size() - 1; i >= 0; --i) {
		VolumedFile f(degree_sorted[i]);
		InputBuffer<Edge> data(f);
		*message_stream << "Computing vertex cover bucket " << i + 1 << "/" << degree_sorted.size() << " edges=" << data.size()
			<< " key_begin=" << degree_sorted[i].key_begin() << " key_end=" << degree_sorted[i].key_end() << endl;
		*message_stream << "Queue nodes=" << queue.size() << " edges=" << edges_queued << " top degree=" << (queue.empty() ? 0 : queue.top().degree) << endl;
		f.remove();
		const OId next_degree = i > 0 ? degree_sorted[i - 1].key_end() - 1 : 0;
#ifdef NDEBUG
		ips4o::parallel::sort(data.begin(), data.end());
#else
		std::sort(data.begin(), data.end());
#endif
		auto it = merge_keys(data.begin(), data.end(), Edge::GetKey());
		while (it.good()) {
			if (cfg.clustering->at(it.key()) != numeric_limits<OId>::max()) {
				++it;
				continue;
			}
			PotentialRep r(it.key());
			r.members.reserve(it.count());
			for (auto j = it.begin(); j != it.end(); ++j) {
				const bool unassigned = cfg.clustering->at(j->node2) == numeric_limits<OId>::max();
				if (!config.no_gvc_reassign || unassigned) {
					r.members.emplace_back(j->node2, j->weight);
					if (unassigned)
						++r.degree;
				}
			}
			edges_queued += r.members.size();
			queue.push(std::move(r));
			++it;
		}
		greedy_vertex_cover(*cfg.clustering, weights, queue, edges_queued, next_degree);
	}
	rmdir(cfg.tmp_dir + "degree_sorted");
	*message_stream << "Time (computing vertex cover): " << timer.seconds() << 's' << endl;
}

void greedy_vertex_cover(Cfg& cfg) {
	const double cov = std::max(config.query_or_target_cover, config.member_cover.get(Cluster::DEFAULT_MEMBER_COVER));
	const bool triplets = config.edge_format == "triplet", symmetric = config.symmetric;
	const bool numeric_ids = config.max_oid.present();
	if (!triplets && symmetric)
		throw runtime_error("--symmetric requires triplet edge format");
	*message_stream << "Coverage cutoff: " << cov << '%' << endl;
	*message_stream << "Edge format: " << (triplets ? "triplet" : "quintuplet") << endl;
	*message_stream << "Symmetric: " << (symmetric ? "yes" : "no") << endl;
	TaskTimer timer;
	unordered_map<Acc, OId> acc2oid;
	OId max_oid;
	size_t db_size;
	if (numeric_ids) {
		if (config.max_oid < 0)
			throw runtime_error("--max-oid must be nonnegative");
		max_oid = config.max_oid;
		if (max_oid >= numeric_limits<size_t>::max())
			throw runtime_error("--max-id exceeds supported maximum");
		db_size = static_cast<size_t>(max_oid) + 1;
		*message_stream << "#Numeric vertex OIDs: " << db_size << endl;
	}
	else {
		config.database.require();
		timer.go("Reading mapping file");
		acc2oid.reserve(File::count_lines(config.database));
		File mapping_file(config.database, "rb");
		string query;
		const char* line;
		while (line = mapping_file.getline(), line[0] != '\0' || !mapping_file.eof()) {
			Util::String::Tokenizer<Util::String::CharDelimiter>(line, Util::String::CharDelimiter('\t')) >> query;
			auto e = acc2oid.emplace(query, acc2oid.size());
			if (!e.second)
				throw runtime_error("Duplicate sequence id found in database file");
		}
		mapping_file.close();
		timer.finish();
		*message_stream << "#Sequences in database: " << acc2oid.size() << endl;
		if (acc2oid.empty())
			throw runtime_error("Database file is empty");
		if (acc2oid.size() > (size_t)numeric_limits<OId>::max())
			throw runtime_error("Input count exceeds supported maximum.");
		db_size = acc2oid.size();
		max_oid = db_size - 1;
	}

	if (cfg.tmp_dir.empty())
		cfg.tmp_dir = config.tmpdir = create_temp_directory(config.tmpdir, "diamond-tmp-") + PATH_SEPARATOR;
	const string base_dir = cfg.tmp_dir;
	//mkdir(base_dir);
	const EdgeInput input{ max_oid, triplets, symmetric, cov, numeric_ids ? nullptr : &acc2oid, cfg.edges };
	RadixedTable rep_sorted = edge_pass_one(base_dir + "rep_sorted" + PATH_SEPARATOR, input, edge_pass_zero(input));
	const DegreePartition p = edge_pass_two(rep_sorted);
	RadixedTable degree_sorted = edge_pass_three(rep_sorted, p, cfg);
	edge_pass_four(degree_sorted, db_size, cfg);
	rmdir(cfg.tmp_dir);

	if (merge_recursive) {
		timer.go("Computing transitive closure");
		for (OId i = 0; i < cfg.clustering->size();) {
			if (cfg.clustering->at(i) != numeric_limits<OId>::max() && cfg.clustering->at(cfg.clustering->at(i)) != cfg.clustering->at(i))
				cfg.clustering->at(i) = cfg.clustering->at(cfg.clustering->at(i));
			else
				++i;
		}
		timer.finish();
	}

	for (size_t i = 0; i < cfg.clustering->size(); ++i)
		if (cfg.clustering->at(i) == numeric_limits<OId>::max())
			cfg.clustering->at(i) = (OId)i;

	std::pmr::monotonic_buffer_resource pool;
	std::pmr::vector<std::pmr::string> acc(&pool);
	if (!numeric_ids) {
		timer.go("Building reverse mapping");
		acc.resize(acc2oid.size());
		for (const auto& i : acc2oid)
			acc.at(i.second) = i.first;
	}
	acc2oid.clear();

	if (!config.centroid_out.empty() || !config.output_file.empty()) {
		timer.go("Generating output");
		OId reps = 0;
		unique_ptr<ofstream> centroid_out;
		if (!config.centroid_out.empty())
			centroid_out.reset(new ofstream(config.centroid_out));
		unique_ptr<ofstream> out;
		if (!config.output_file.empty())
			out.reset(new ofstream(config.output_file));
		for (size_t i = 0; i < cfg.clustering->size(); ++i) {
			if (cfg.clustering->at(i) == i) {
				++reps;
				if (!config.centroid_out.empty()) {
					if (numeric_ids)
						*centroid_out << i << endl;
					else
						*centroid_out << acc[i] << endl;
				}
			}
			if (!config.output_file.empty()) {
				if (numeric_ids)
					*out << cfg.clustering->at(i) << '\t' << i << endl;
				else
					*out << acc[cfg.clustering->at(i)] << '\t' << acc[i] << endl;
			}
		}
		if (centroid_out)
			centroid_out->close();
		if (out)
			out->close();
		timer.finish();
		*message_stream << "#Clusters: " << reps << endl;
	}
}

void greedy_vertex_cover() {
	Cfg cfg;
	greedy_vertex_cover(cfg);
}

}