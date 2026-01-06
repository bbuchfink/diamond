/****
DIAMOND protein aligner
Copyright (C) 2021-2024 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#include <fstream>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <memory>
#include "util/tsv/tsv.h"
#include "basic/config.h"
#include "util/log_stream.h"
#include "util/string/tokenizer.h"
#include "util/algo/algo.h"
#include "util/system/system.h"
#include "cluster/cluster.h"

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
//using Acc = FixedString<32>;
using Acc = string;

void greedy_vertex_cover() {
	config.database.require();
	using Int = uint64_t;
	using Edge = Util::Algo::Edge<Int>;
	const double cov = std::max(config.query_or_target_cover, config.member_cover.get(Cluster::DEFAULT_MEMBER_COVER));
	const bool triplets = config.edge_format == "triplet", symmetric = config.symmetric;
	if (!triplets && symmetric)
		throw std::runtime_error("--symmetric requires triplet edge format");
	message_stream << "Coverage cutoff: " << cov << '%' << endl;
	TaskTimer timer("Reading mapping file");
	//unordered_map<Acc, OId, Acc::Hash> acc2oid;
	unordered_map<Acc, OId> acc2oid;
	acc2oid.reserve(Util::Tsv::count_lines(config.database));
	TextInputFile mapping_file(config.database);
	string query;
	while (mapping_file.getline(), !mapping_file.line.empty() || !mapping_file.eof()) {
		Util::String::Tokenizer<Util::String::CharDelimiter>(mapping_file.line, Util::String::CharDelimiter('\t')) >> query;
		acc2oid.emplace(query, acc2oid.size());
	}
	mapping_file.close();
	timer.finish();
	message_stream << "#OIds: " << acc2oid.size() << endl;
	if (acc2oid.size() > (size_t)numeric_limits<Int>::max())
		throw runtime_error("Input count exceeds supported maximum.");

	timer.go("Counting input lines");
	atomic<int64_t> lines(0);
	function<void(int64_t, const char*, const char*)> fn([&](int64_t, const char* begin, const char* end) {
		Util::String::LineIterator it(begin, end);
		int64_t n = 0;
		double qcov, tcov;
		while (it.good()) {
			string line = *it;
			if (triplets)
				if (symmetric)
					n += 2;
				else
					++n;
			else {
				Util::String::Tokenizer<Util::String::CharDelimiter>(line, Util::String::CharDelimiter('\t')) >> Util::String::Skip() >> Util::String::Skip() >> qcov >> tcov;
				if (qcov >= cov)
					++n;
				if (tcov >= cov)
					++n;
			}
			++it;
		}
		lines += n;
		});
	Util::Tsv::File(Util::Tsv::Schema(), config.edges).read(INT64_MAX, config.threads_, fn);
	timer.finish();
	message_stream << "#Lines: " << lines << endl;

	timer.go("Allocating memory");
	vector<Edge> edges;
	mutex mtx;
	edges.reserve(lines);

	timer.go("Reading input lines");
	function<void(int64_t, const char*, const char*)> fn2([&](int64_t, const char* begin, const char* end) {
		Util::String::LineIterator it(begin, end);
		vector<Edge> e;
		string query, target;
		double qcov, tcov, evalue;
		while (it.good()) {
			string line = *it;
			Util::String::Tokenizer<Util::String::CharDelimiter> tok(line, Util::String::CharDelimiter('\t'));
			tok >> query >> target;
			if (!triplets)
				tok >> qcov >> tcov;
			tok >> evalue;
			if (triplets || tcov >= cov || qcov >= cov) {
				const auto q = acc2oid.at(query), t = acc2oid.at(target);
				if (q == t) {
					++it;
					continue;
				}
				if (triplets) {
					e.emplace_back(t, q, evalue);
					if (symmetric)
						e.emplace_back(q, t, evalue);
				}
				else {
					if (tcov >= cov)
						e.emplace_back(q, t, evalue);
					if (qcov >= cov)
						e.emplace_back(t, q, evalue);
				}
			}
			++it;
		}
		{
			lock_guard<mutex> lock(mtx);
			edges.insert(edges.end(), e.cbegin(), e.cend());
		}
		});
	Util::Tsv::File(Util::Tsv::Schema(), config.edges).read(INT64_MAX, config.threads_, fn2);
	timer.finish();
	log_rss();

	timer.go("Making flat array");
	FlatArray<Edge> edge_array = make_flat_array_dense(std::move(edges), (Int)acc2oid.size(), config.threads_, Edge::GetKey());
	timer.finish();
	log_rss();

	auto r = Util::Algo::greedy_vertex_cover<Int>(edge_array, nullptr, !config.strict_gvc, !config.no_gvc_reassign, (Int)atoi(config.connected_component_depth.front().c_str()));

	timer.go("Building reverse mapping");
	vector<string> acc(acc2oid.size());
	for (const auto& i : acc2oid)
		//acc[i.second] = string(i.first.chars.data());
		acc[i.second] = i.first;

	timer.go("Generating output");
	int64_t c = 0;
	unique_ptr<ofstream> centroid_out;
	if (!config.centroid_out.empty())
		centroid_out.reset(new ofstream(config.centroid_out));
	unique_ptr<ofstream> out;
	if (!config.output_file.empty())
		out.reset(new ofstream(config.output_file));
	for (int64_t i = 0; i < (int64_t)r.size(); ++i) {
		if (r[i] == i) {
			++c;
			if (!config.centroid_out.empty())
				*centroid_out << acc[i] << endl;
		}
		if (!config.output_file.empty())
			*out << acc[r[i]] << '\t' << acc[i] << endl;
	}
	if(centroid_out)
		centroid_out->close();
	if (out)
		out->close();
	timer.finish();
	message_stream << "#Centroids: " << c << endl;
}