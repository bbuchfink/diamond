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
#include "util/algo/algo.h"
#include "util/system/system.h"
#include "cluster/cluster.h"
#include "cluster/multinode/file_array.h"
#include "cluster/multinode/input_buffer.h"
#include "util/memory/memory_resource.h"
#include "tools.h"
#include "data/fasta/parser.h"

namespace Util { namespace Algo {
inline void serialize(const Edge<uint64_t>& e, CompressedBuffer& buf) {
	buf.write(e.node1);
	buf.write(e.node2);
	buf.write(e.weight);
}
}}

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
using Int = uint64_t;
using Edge = Util::Algo::Edge<Int>;

namespace GVC {

static const bool merge_recursive = true;
static const uint64_t RADIX_COUNT = 256;
static const int RADIX_BITS = 8;

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

static RadixedTable edge_pass_one(const string& base_dir, const OId max_oid, bool triplets, bool symmetric, double cov, const unordered_map<Acc, OId>& acc2oid) {
	mkdir(base_dir);
	FileArray file_array(base_dir, RADIX_COUNT, 0, true);
	const int shift = std::max(bit_length(max_oid) - RADIX_BITS, 0);
	atomic<uint64_t> total_lines(0);

	TaskTimer timer("Reading input lines");
	function<void(int64_t, const char*, const char*)> fn2([&](int64_t, const char* begin, const char* end) {
		size_t line_count = 0;
		Util::String::LineIterator it(begin, end);
		thread_local BufferArray buffers(file_array, RADIX_COUNT);
		string query, target;
		double qcov, tcov, weight;
		auto emit = [&](Int n1, Int n2, double w) {
			Edge edge(n1, n2, w);
			const uint64_t radix = (uint64_t)n1 >> shift;
			buffers.write(radix, &edge, 1, 1);
			};
		while (it.good()) {
			++line_count;
			string line = *it;
			if (ends_with(line, "\r"))
				line.pop_back();
			Util::String::Tokenizer<Util::String::CharDelimiter> tok(line, Util::String::CharDelimiter('\t'));
			tok >> query >> target;
			if (!triplets)
				tok >> qcov >> tcov;
			tok >> weight;
			if (triplets || tcov >= cov || qcov >= cov) {
				const auto q = acc2oid.at(query), t = acc2oid.at(target);
				if (q == t) {
					++it;
					continue;
				}
				if (triplets) {
					emit(t, q, weight);
					if (symmetric)
						emit(q, t, weight);
				}
				else {
					if (tcov >= cov)
						emit(q, t, weight);
					if (qcov >= cov)
						emit(t, q, weight);
				}
			}
			++it;
		}
		total_lines.fetch_add(line_count, std::memory_order_relaxed);
		});
	File in(config.edges, "rb");
	in.read_text_mt(INT64_MAX, config.threads_, fn2);
	timer.finish();
	log_rss();
	*message_stream << "#Input lines: " << total_lines.load(std::memory_order_relaxed) << endl;
	*message_stream << "#Input edges: " << file_array.records_total() << endl;
	file_array.close();
	return file_array.buckets(shift);
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
	unordered_map<OId, OId> degrees;
	for (RadixedTable::const_iterator it = rep_sorted.begin(); it != rep_sorted.end(); ++it) {
		VolumedFile f(*it);
		InputBuffer<Edge> data(f);
		*message_stream << "Finding neighbor counts bucket " << it - rep_sorted.begin() + 1 << "/" << rep_sorted.size() << " records=" << data.size() << endl;
		ips4o::parallel::sort(data.begin(), data.end());
		auto i = merge_keys(data.begin(), data.end(), Edge::GetKey());
		while (i.good()) {
			const OId d = node_degree(i.begin(), i.end());
			degrees[d] += d;
			++i;
		}
	}
	return DegreePartition(degrees, RADIX_COUNT);
}

static RadixedTable edge_pass_three(const RadixedTable& rep_sorted, const DegreePartition& p, Cfg& cfg) {
	const string base_dir = cfg.tmp_dir + "degree_sorted";
	mkdir(base_dir);
	FileArray file_array(base_dir, p.size(), 0, true);
	BufferArray buffers(file_array, p.size());
	for (RadixedTable::const_iterator it = rep_sorted.begin(); it != rep_sorted.end(); ++it) {
		VolumedFile f(*it);
		InputBuffer<Edge> data(f);
		*message_stream << "Writing degree sorted edges " << it - rep_sorted.begin() + 1 << "/" << rep_sorted.size() << " records=" << data.size() << endl;
		ips4o::parallel::sort(data.begin(), data.end());
		auto i = merge_keys(data.begin(), data.end(), Edge::GetKey());
		while (i.good()) {
			const int bucket = p.bucket_index(node_degree(i.begin(), i.end()));
			for (auto j = i.begin(); j != i.end(); ++j) {
				if (j->node2 != j->node1 && (j == i.begin() || j->node2 != (j - 1)->node2))
					buffers.write(bucket, *j);
			}
			++i;
		}
		f.remove();
	}
	buffers.finish();
	file_array.close();
	rmdir(cfg.tmp_dir + "rep_sorted");
	return file_array.buckets(p);
}

static vector<OId> edge_pass_four(const RadixedTable& degree_sorted, OId db_size, Cfg& cfg) {
	vector<OId> clustering(db_size, numeric_limits<OId>::max());
	vector<double> weights(db_size);
	RepQueue queue;
	uint64_t edges_queued = 0;
	for (int i = degree_sorted.size()  - 1; i >= 0; --i) {
		VolumedFile f(degree_sorted[i]);
		InputBuffer<Edge> data(f);
		*message_stream << "Computing vertex cover bucket " << i + 1 << "/" << degree_sorted.size() << " edges=" << data.size()
			<< " key_begin=" << degree_sorted[i].key_begin() << " key_end=" << degree_sorted[i].key_end() << endl;
		*message_stream << "Queue nodes=" << queue.size() << " edges=" << edges_queued << " top degree=" << (queue.empty() ? 0 : queue.top().degree) << endl;
		f.remove();
		const OId next_degree = i > 0 ? degree_sorted[i - 1].key_end() - 1 : 0;
		ips4o::parallel::sort(data.begin(), data.end());
		auto it = merge_keys(data.begin(), data.end(), Edge::GetKey());
		while (it.good()) {
			if (clustering[it.key()] != numeric_limits<OId>::max()) {
				++it;
				continue;
			}
			PotentialRep r(it.key());
			r.members.reserve(it.count());
			for (auto j = it.begin(); j != it.end(); ++j) {
				const bool unassigned = clustering[j->node2] == numeric_limits<OId>::max();
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
		greedy_vertex_cover(clustering, weights, queue, edges_queued, next_degree);
	}
	rmdir(cfg.tmp_dir + "degree_sorted");
	return clustering;
}

void greedy_vertex_cover(Cfg& cfg) {
	config.database.require();
	const double cov = std::max(config.query_or_target_cover, config.member_cover.get(Cluster::DEFAULT_MEMBER_COVER));
	const bool triplets = config.edge_format == "triplet", symmetric = config.symmetric;
	if (!triplets && symmetric)
		throw runtime_error("--symmetric requires triplet edge format");
	*message_stream << "Coverage cutoff: " << cov << '%' << endl;
	*message_stream << "Edge format: " << (triplets ? "triplet" : "quintuplet") << endl;
	*message_stream << "Symmetric: " << (symmetric ? "yes" : "no") << endl;
	TaskTimer timer("Reading mapping file");
	unordered_map<Acc, OId> acc2oid;
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
	if (acc2oid.size() > (size_t)numeric_limits<Int>::max())
		throw runtime_error("Input count exceeds supported maximum.");
	const OId max_oid = acc2oid.size() - 1;

	if (cfg.tmp_dir.empty())
		cfg.tmp_dir = config.tmpdir = create_temp_directory(config.tmpdir, "diamond-tmp-") + PATH_SEPARATOR;
	const string base_dir = cfg.tmp_dir;
	//mkdir(base_dir);
	RadixedTable rep_sorted = edge_pass_one(base_dir + "rep_sorted" + PATH_SEPARATOR, max_oid, triplets, symmetric, cov, acc2oid);
	const DegreePartition p = edge_pass_two(rep_sorted);
	RadixedTable degree_sorted = edge_pass_three(rep_sorted, p, cfg);
	vector<OId> clustering = edge_pass_four(degree_sorted, acc2oid.size(), cfg);
	rmdir(cfg.tmp_dir);

	if (merge_recursive) {
		timer.go("Computing transitive closure");
		for (OId i = 0; i < clustering.size();) {
			if (clustering[i] != numeric_limits<OId>::max() && clustering[clustering[i]] != clustering[i])
				clustering[i] = clustering[clustering[i]];
			else
				++i;
		}
	}

	timer.go("Building reverse mapping");
	std::pmr::monotonic_buffer_resource pool;
	std::pmr::vector<std::pmr::string> acc(acc2oid.size(), &pool);
	for (const auto& i : acc2oid)
		acc.at(i.second) = i.first;
	acc2oid.clear();

	timer.go("Generating output");
	OId reps = 0;
	unique_ptr<ofstream> centroid_out;
	if (!config.centroid_out.empty())
		centroid_out.reset(new ofstream(config.centroid_out));
	unique_ptr<ofstream> out;
	if (!config.output_file.empty())
		out.reset(new ofstream(config.output_file));
	for (size_t i = 0; i < clustering.size(); ++i) {
		if (clustering[i] == numeric_limits<OId>::max())
			clustering[i] = i;
		if (clustering[i] == i) {
			++reps;
			if (!config.centroid_out.empty())
				*centroid_out << acc[i] << endl;
		}
		if (!config.output_file.empty())
			*out << acc[clustering[i]] << '\t' << acc[i] << endl;
	}
	if (centroid_out)
		centroid_out->close();
	if (out)
		out->close();
	timer.finish();
	*message_stream << "#Clusters: " << reps << endl;
}

void greedy_vertex_cover() {
	Cfg cfg;
	greedy_vertex_cover(cfg);
}

}