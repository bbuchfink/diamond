/****
Copyright © 2013-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-2-Clause

#include <sstream>
#include "cluster.h"
#include "util/tsv/tsv.h"
#include "util/string/tokenizer.h"
#include "basic/config.h"
#include "util/log_stream.h"

const char* const HEADER_LINE = "centroid\tmember";

using std::pair;
using std::endl;
using std::vector;
using std::ofstream;
using std::ostream;
using std::runtime_error;
using std::string;
using std::for_each;
using std::floor;
using std::unique_ptr;
using std::back_inserter;
using std::less;
using std::stringstream;
using namespace Util::Tsv;

namespace Cluster {

template<typename Int>
pair<FlatArray<Int>, vector<Int>> read(const string& file_name, const SequenceFile& db, CentroidSorted) {
	const int64_t lines = Util::Tsv::count_lines(file_name);
	TextInputFile in(file_name);
	string centroid, member;
	vector<pair<Int, Int>> pairs;
	pairs.reserve(lines);
	int64_t mappings = 0;
	if (TabularFormat::header_format(::Config::cluster) == Header::SIMPLE) {
		in.getline();
		if (in.line != HEADER_LINE)
			throw runtime_error("Clusters file is missing header line.");
	}
	while (in.getline(), !in.eof() || !in.line.empty()) {
		Util::String::Tokenizer<Util::String::CharDelimiter>(in.line, Util::String::CharDelimiter('\t')) >> centroid >> member;
		const Int centroid_oid = (Int)db.accession_to_oid(centroid).front();
		const Int member_oid = (Int)db.accession_to_oid(member).front();
		pairs.emplace_back(centroid_oid, member_oid);
		++mappings;
		if (mappings % 1000000 == 0)
			log_stream << "#Entries: " << mappings << endl;
	}
	in.close();
	return make_flat_array(pairs.begin(), pairs.end(), config.threads_);
}

template pair<FlatArray<uint32_t>, vector<uint32_t>> read(const string&, const SequenceFile&, CentroidSorted);
template pair<FlatArray<uint64_t>, vector<uint64_t>> read(const string&, const SequenceFile&, CentroidSorted);

template<typename Int>
vector<Int> read(const string& file_name, const SequenceFile& db) {
	const int64_t lines = Util::Tsv::count_lines(file_name);
	TextInputFile in(file_name);
	string centroid, member;
	vector<Int> v(db.sequence_count());
	uint64_t mappings = 0;
	if (TabularFormat::header_format(::Config::cluster) == Header::SIMPLE) {
		in.getline();
		if (in.line != HEADER_LINE)
			throw runtime_error("Clustering input file is missing header line.");
	}
	while (in.getline(), !in.eof() || !in.line.empty()) {
		Util::String::Tokenizer<Util::String::CharDelimiter>(in.line, Util::String::CharDelimiter('\t')) >> centroid >> member;
		const auto centroid_oid = db.accession_to_oid(centroid);
		const auto member_oid = db.accession_to_oid(member);
		v[member_oid.front()] = (Int)centroid_oid.front();
		++mappings;
		if (mappings % 1000000 == 0)
			log_stream << "#Entries: " << mappings << endl;
	}
	in.close();
	if (mappings != db.sequence_count())
		throw runtime_error("Invalid/incomplete clustering.");
	return v;
}

template vector<uint32_t> read(const string&, const SequenceFile&);
template vector<uint64_t> read(const string&, const SequenceFile&);

template<typename Int>
vector<Int> member2centroid_mapping(const FlatArray<Int>& clusters, const vector<Int>& centroids) {
	vector<Int> v(clusters.data_size());
	for (Int i = 0; i < (Int)centroids.size(); ++i) {
		for (Int j : clusters[i])
			v[j] = centroids[i];
	}
	return v;
}
template vector<uint64_t> member2centroid_mapping(const FlatArray<uint64_t>&, const vector<uint64_t>&);

template<typename Int>
pair<FlatArray<Int>, vector<Int>> cluster_sorted(const vector<Int>& mapping) {
	vector<pair<Int, Int>> pairs;
	pairs.reserve(mapping.size());
	for (typename vector<Int>::const_iterator i = mapping.begin(); i < mapping.end(); ++i)
		pairs.emplace_back(*i, Int(i - mapping.begin()));
	return make_flat_array(pairs.begin(), pairs.end(), config.threads_);
}

template pair<FlatArray<uint32_t>, vector<uint32_t>> cluster_sorted<uint32_t>(const vector<uint32_t>&);
template pair<FlatArray<uint64_t>, vector<uint64_t>> cluster_sorted<uint64_t>(const vector<uint64_t>&);

void output(File& out, SequenceFile& db, File& oid_to_centroid_oid) {
	unique_ptr<File> sorted1(oid_to_centroid_oid.sort(1, config.threads_));
	unique_ptr<File> joined1(join(*sorted1, db.seqid_file(), 1, 0, { {0,0},{1,1} }));
	sorted1.reset();
	unique_ptr<File> sorted2(joined1->sort(0, config.threads_));
	join(*sorted2, db.seqid_file(), 0, 0, { {1,1}, {0,1} }, out);
}

template<typename Int>
void output_mem(File& out, SequenceFile& db, const FlatArray<Int>& clusters, const vector<Int>& centroids) {
	if (config.oid_output) {
		for (Int i = 0; i < (Int)centroids.size(); ++i) {
			const string centroid = std::to_string(centroids[i]);
			for (auto j = clusters.cbegin(i); j != clusters.cend(i); ++j)
				out.write_record(centroid, std::to_string(*j));
		}
	}
	else {
		const Util::Tsv::Table acc_mapping = db.seqid_file().read(config.threads_);
		for (Int i = 0; i < (Int)centroids.size(); ++i) {
			const string centroid = acc_mapping[centroids[i]].template get<string>(0);
			for (auto j = clusters.cbegin(i); j != clusters.cend(i); ++j)
				out.write_record(centroid, acc_mapping[*j].template get<string>(0));
		}
	}
}

template void output_mem<uint32_t>(File&, SequenceFile&, const FlatArray<uint32_t>&, const vector<uint32_t>&);
template void output_mem<uint64_t>(File&, SequenceFile&, const FlatArray<uint64_t>&, const vector<uint64_t>&);

template<typename Int>
void output_mem(File& out, SequenceFile& db, const vector<Int>& mapping) {
	vector<Int> centroids;
	FlatArray<Int> clusters;
	tie(clusters, centroids) = cluster_sorted(mapping);
	output_mem<Int>(out, db, clusters, centroids);
}

template void output_mem<uint32_t>(File&, SequenceFile&, const vector<uint32_t>&);
template void output_mem<uint64_t>(File&, SequenceFile&, const vector<uint64_t>&);

template<typename Int>
void output_mem(File& out, SequenceFile& db, File& oid_to_centroid_oid) {
	vector<pair<Int, Int>> centroid_oid;
	oid_to_centroid_oid.template read<Int, Int>(back_inserter(centroid_oid));
	vector<Int> centroids;
	FlatArray<Int> clusters;
	tie(clusters, centroids) = make_flat_array(centroid_oid.begin(), centroid_oid.end(), config.threads_);
	output_mem<Int>(out, db, clusters, centroids);
}

void output_mem(File& out, SequenceFile& db, File& oid_to_centroid_oid) {
	if (db.sequence_count() > (int64_t)INT32_MAX)
		output_mem<uint64_t>(out, db, oid_to_centroid_oid);
	else
		output_mem<uint32_t>(out, db, oid_to_centroid_oid);
}

template<typename Int>
pair<vector<Int>, vector<Int>> split(const vector<Int>& mapping) {
	pair<vector<Int>, vector<Int>> r;
	r.second.reserve(mapping.size());
	for (Int i = 0; i < (Int)mapping.size(); ++i)
		if (mapping[i] == i)
			r.first.push_back(i);
		else
			r.second.push_back(i);
	return r;
}

template pair<vector<uint32_t>, vector<uint32_t>> split(const vector<uint32_t>&);
template pair<vector<uint64_t>, vector<uint64_t>> split(const vector<uint64_t>&);

vector<SuperBlockId> member_counts(const vector<SuperBlockId>& mapping) {
	vector<SuperBlockId> v(mapping.size(), 0);
	for (SuperBlockId c : mapping)
		++v[c];
	return v;
}

void init_thresholds() {
	if (config.member_cover.present() && config.mutual_cover.present())
		throw runtime_error("--member-cover and --mutual-cover are mutually exclusive.");
	if (!config.mutual_cover.present())
		config.member_cover.set_if_blank(DEFAULT_MEMBER_COVER);
	if (!config.approx_min_id.present())
		config.approx_min_id = config.command == ::Config::DEEPCLUST ? 0.0 : (config.command == ::Config::LINCLUST ? 90.0 : 50.0);
	if (config.soft_masking.empty())
		config.soft_masking = "tantan";
	if (!config.masking_.present())
		config.masking_ = "0";
	if (config.approx_min_id < 90.0 || config.mutual_cover.present())
		return;
	config.diag_filter_id.set_if_blank(config.approx_min_id - (config.command == ::Config::CLUSTER_REASSIGN ? 10.0 : 10.0));
	if (config.approx_min_id < 90.0)
		return;
	if(config.command == ::Config::CLUSTER_REASSIGN) {
		config.diag_filter_cov.set_if_blank(config.member_cover > 50.0 ? config.member_cover - 10.0 : 0.0);
	}
	else {
		config.diag_filter_cov.set_if_blank(config.member_cover > 50.0 ? config.member_cover - 10.0 : 0.0);
	}
}

File* open_out_tsv() {
	File* file = new File(Schema{ Type::STRING, Type::STRING }, config.output_file, Flags::WRITE);
	if (TabularFormat::header_format(::Config::cluster) == Header::SIMPLE)
		file->write_record("centroid", "member");
	return file;
}

vector<BlockId> len_sorted_clust(const FlatArray<Util::Algo::Edge<SuperBlockId>>& edges) {
	static constexpr BlockId NIL = std::numeric_limits<BlockId>::max();
	vector<BlockId> v(edges.size(), NIL);
	for (uint64_t i = 0; i < edges.size(); ++i) {
		if (v[i] != NIL)
			continue;
		v[i] = (BlockId)i;
		for (auto it = edges.cbegin(i); it != edges.cend(i); ++it)
			if (v[it->node2] == NIL)
				v[it->node2] = (BlockId)i;
	}
	return v;
}

void output_edges(const string& file, SequenceFile& db, const vector<Util::Algo::Edge<SuperBlockId>>& edges) {
	File out(Schema{ Type::STRING, Type::STRING }, file, Flags::WRITE);
	const Util::Tsv::Table acc_mapping = db.seqid_file().read(config.threads_);
	for (const auto& e : edges) {
		out.write_record(acc_mapping[e.node1].get<string>(0), acc_mapping[e.node2].get<string>(0));
	}
}

double round_value(const vector<string>& par, const string& name, int round, int round_count) {
	if (par.empty())
		return 0.0;
	if (round >= round_count - 1)
		return 0.0;
	if ((int64_t)par.size() >= round_count)
		throw runtime_error("Too many values provided for " + name);
	vector<double> v;
	for(const string& s : par) {
		stringstream ss(s);
		double i;
		ss >> i;
		if (ss.fail() || !ss.eof())
			throw runtime_error("Invalid value provided for " + name + ": " + s);
		v.push_back(i);
	}
	v.insert(v.begin(), round_count - 1 - v.size(), v.front());
	return v[round];
}

}