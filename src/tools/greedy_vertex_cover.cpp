#include <iostream>
#include <fstream>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include "../util/tsv/tsv.h"
#include "../util/tsv/file.h"
#include "../basic/config.h"
#include "../util/log_stream.h"
#include "../util/string/fixed_string.h"
#include "../util/string/tokenizer.h"
#include "../util/algo/algo.h"
#include "../util/system/system.h"

using std::ofstream;
using std::mutex;
using std::lock_guard;
using std::endl;
using std::string;
using std::atomic;
using std::vector;
using std::unordered_map;
using std::cout;
using std::numeric_limits;
using std::runtime_error;
//using Acc = FixedString<32>;
using Acc = string;

void greedy_vertex_cover() {
	using Int = int64_t;
	using Edge = Util::Algo::Edge<Int>;
	const double cov = std::max(config.query_or_target_cover, config.member_cover);
	const bool triplets = config.edge_format == "triplet";
	message_stream << "Coverage cutoff: " << cov << '%' << endl;
	task_timer timer("Reading mapping file");
	//unordered_map<Acc, OId, Acc::Hash> acc2oid;
	unordered_map<Acc, OId> acc2oid;
	acc2oid.reserve(Util::Tsv::count_lines(config.database));
	TextInputFile mapping_file(config.database);
	string query;
	while (mapping_file.getline(), !mapping_file.line.empty() || !mapping_file.eof()) {
		Util::String::Tokenizer(mapping_file.line, "\t") >> query;
		acc2oid.emplace(query, acc2oid.size());
	}
	mapping_file.close();
	timer.finish();
	message_stream << "#OIds: " << acc2oid.size() << endl;
	if (acc2oid.size() > (size_t)numeric_limits<Int>::max())
		throw runtime_error("Input count exceeds supported maximum.");

	timer.go("Counting input lines");
	atomic<int64_t> lines(0);
	std::function<void(int64_t, const char*, const char*)> fn([&](int64_t, const char* begin, const char* end) {
		Util::Tsv::LineIterator it(begin, end);
		int64_t n = 0;
		double qcov, tcov;
		while (it.good()) {
			string line = *it;
			if (triplets)
				++n;
			else {
				Util::String::Tokenizer(line, "\t") >> Util::String::Skip() >> Util::String::Skip() >> qcov >> tcov;
				if (qcov >= cov)
					++n;
				if (tcov >= cov)
					++n;
			}
			++it;
		}
		lines += n;
		});
	Util::Tsv::File(Util::Tsv::Schema(), config.edges).read(config.threads_, fn);
	timer.finish();
	message_stream << "#Lines: " << lines << endl;

	timer.go("Allocating memory");
	vector<Edge> edges;
	mutex mtx;
	edges.reserve(lines);

	timer.go("Reading input lines");
	std::function<void(int64_t, const char*, const char*)> fn2([&](int64_t, const char* begin, const char* end) {
		Util::Tsv::LineIterator it(begin, end);
		vector<Edge> e;
		string query, target;
		double qcov, tcov, evalue;
		while (it.good()) {
			string line = *it;
			Util::String::Tokenizer tok(line, "\t");
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
				if(triplets)
					e.emplace_back(t, q, evalue);
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
	Util::Tsv::File(Util::Tsv::Schema(), config.edges).read(config.threads_, fn2);
	timer.finish();
	log_rss();

	timer.go("Making flat array");
	FlatArray<Edge> edge_array = make_flat_array_dense(move(edges), (Int)acc2oid.size(), config.threads_, Edge::GetKey());
	timer.finish();
	log_rss();

	auto r = Util::Algo::greedy_vertex_cover(edge_array, nullptr, !config.strict_gvc);

	timer.go("Building reverse mapping");
	vector<string> acc(acc2oid.size());
	for (const auto& i : acc2oid)
		//acc[i.second] = string(i.first.chars.data());
		acc[i.second] = i.first;

	timer.go("Generating output");
	int64_t c = 0;
	ofstream centroid_out;
	if (!config.centroid_out.empty())
		centroid_out = ofstream(config.centroid_out);
	ofstream out;
	if (!config.output_file.empty())
		out = ofstream(config.output_file);
	for (int64_t i = 0; i < (int64_t)r.size(); ++i) {
		if (r[i] == i) {
			++c;
			if (!config.centroid_out.empty())
				centroid_out << acc[i] << endl;
		}
		if (!config.output_file.empty())
			out << acc[r[i]] << '\t' << acc[i] << endl;
	}
	centroid_out.close();
	out.close();
	timer.finish();
	message_stream << "#Centroids: " << c << endl;
}