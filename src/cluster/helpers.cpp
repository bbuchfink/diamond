/****
DIAMOND protein aligner
Copyright (C) 2021-2023 Max Planck Society for the Advancement of Science e.V.

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
#include "cluster.h"
#include "../util/tsv/tsv.h"
#include "../util/string/tokenizer.h"
#include "../basic/config.h"
#include "../search/search.h"
#define _REENTRANT
#include "../lib/ips4o/ips4o.hpp"

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
	if (Blast_tab_format::header_format(::Config::cluster) == Header::SIMPLE) {
		in.getline();
		if (in.line != HEADER_LINE)
			throw runtime_error("Clusters file is missing header line.");
	}
	while (in.getline(), !in.eof() || !in.line.empty()) {
		Util::String::Tokenizer(in.line, "\t") >> centroid >> member;
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

template pair<FlatArray<int32_t>, vector<int32_t>> read(const string&, const SequenceFile&, CentroidSorted);
template pair<FlatArray<int64_t>, vector<int64_t>> read(const string&, const SequenceFile&, CentroidSorted);

template<typename Int>
vector<Int> read(const string& file_name, const SequenceFile& db) {
	const int64_t lines = Util::Tsv::count_lines(file_name);
	TextInputFile in(file_name);
	string centroid, member;
	vector<Int> v(db.sequence_count());
	int64_t mappings = 0;
	if (Blast_tab_format::header_format(::Config::cluster) == Header::SIMPLE) {
		in.getline();
		if (in.line != HEADER_LINE)
			throw runtime_error("Clustering input file is missing header line.");
	}
	while (in.getline(), !in.eof() || !in.line.empty()) {
		Util::String::Tokenizer(in.line, "\t") >> centroid >> member;
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

template vector<int32_t> read(const string&, const SequenceFile&);
template vector<int64_t> read(const string&, const SequenceFile&);

template<typename Int>
vector<Int> member2centroid_mapping(const FlatArray<Int>& clusters, const vector<Int>& centroids) {
	vector<Int> v(clusters.data_size());
	for (Int i = 0; i < (Int)centroids.size(); ++i) {
		for (Int j : clusters[i])
			v[j] = centroids[i];
	}
	return v;
}

template<typename Int>
pair<FlatArray<Int>, vector<Int>> cluster_sorted(const vector<Int>& mapping) {
	vector<pair<Int, Int>> pairs;
	pairs.reserve(mapping.size());
	for (typename vector<Int>::const_iterator i = mapping.begin(); i < mapping.end(); ++i)
		pairs.emplace_back(*i, Int(i - mapping.begin()));
	return make_flat_array(pairs.begin(), pairs.end(), config.threads_);
}

template pair<FlatArray<int32_t>, vector<int32_t>> cluster_sorted<int32_t>(const vector<int32_t>&);
template pair<FlatArray<int64_t>, vector<int64_t>> cluster_sorted<int64_t>(const vector<int64_t>&);

void output(File& out, SequenceFile& db, File& oid_to_centroid_oid) {
	unique_ptr<File> sorted1(oid_to_centroid_oid.sort(1, config.threads_));
	unique_ptr<File> joined1(join(*sorted1, db.seqid_file(), 1, 0, { {0,0},{1,1} }));
	sorted1.reset();
	unique_ptr<File> sorted2(joined1->sort(0, config.threads_));
	join(*sorted2, db.seqid_file(), 0, 0, { {1,1}, {0,1} }, out);
}

template<typename Int>
void output_mem(File& out, SequenceFile& db, const FlatArray<Int>& clusters, const vector<Int>& centroids) {
	const Util::Tsv::Table acc_mapping = db.seqid_file().read(config.threads_);
	for (Int i = 0; i < (Int)centroids.size(); ++i) {
		const string centroid = acc_mapping[centroids[i]].template get<string>(0);
		for(auto j = clusters.cbegin(i); j != clusters.cend(i); ++j)
			out.write_record(centroid, acc_mapping[*j].template get<string>(0));
	}
}

template void output_mem<int32_t>(File&, SequenceFile&, const FlatArray<int32_t>&, const vector<int32_t>&);
template void output_mem<int64_t>(File&, SequenceFile&, const FlatArray<int64_t>&, const vector<int64_t>&);

template<typename Int>
void output_mem(File& out, SequenceFile& db, const vector<Int>& mapping) {
	vector<Int> centroids;
	FlatArray<Int> clusters;
	tie(clusters, centroids) = cluster_sorted(mapping);
	output_mem<Int>(out, db, clusters, centroids);
}

template void output_mem<int32_t>(File&, SequenceFile&, const vector<int32_t>&);
template void output_mem<int64_t>(File&, SequenceFile&, const vector<int64_t>&);

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
		output_mem<int64_t>(out, db, oid_to_centroid_oid);
	else
		output_mem<int32_t>(out, db, oid_to_centroid_oid);
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

template pair<vector<int32_t>, vector<int32_t>> split(const vector<int32_t>&);
template pair<vector<int64_t>, vector<int64_t>> split(const vector<int64_t>&);

vector<SuperBlockId> member_counts(const vector<SuperBlockId>& mapping) {
	vector<SuperBlockId> v(mapping.size(), 0);
	for (SuperBlockId c : mapping)
		++v[c];
	return v;
}

void init_thresholds() {
	if (!config.approx_min_id.present())
		config.approx_min_id = config.command == ::Config::DEEPCLUST ? 0.0 : (config.command == ::Config::LINCLUST ? 90.0 : 50.0);
	if (config.soft_masking.empty())
		config.soft_masking = "tantan";
	if (!config.masking_.present())
		config.masking_ = "0";
	if (config.approx_min_id < 90.0)
		return;
	if(config.command == ::Config::CLUSTER_REASSIGN) {
		config.diag_filter_id = 80.0;
		config.diag_filter_cov = config.member_cover > 50.0 ? config.member_cover - 10.0 : 0.0;
	}
	else {
		config.diag_filter_id = 85.0;
		config.diag_filter_cov = config.member_cover > 50.0 ? config.member_cover - 10.0 : 0.0;
	}
}

File* open_out_tsv() {
	File* file = new File(Schema{ Type::STRING, Type::STRING }, config.output_file, Flags::WRITE);
	if (Blast_tab_format::header_format(::Config::cluster) == Header::SIMPLE)
		file->write_record("centroid", "member");
	return file;
}

vector<BlockId> len_sorted_clust(const FlatArray<Util::Algo::Edge<SuperBlockId>>& edges) {
	vector<BlockId> v(edges.size(), -1);
	for (int64_t i = 0; i < edges.size(); ++i) {
		if (v[i] != -1)
			continue;
		v[i] = (BlockId)i;
		for (auto it = edges.cbegin(i); it != edges.cend(i); ++it)
			if (v[it->node2] == -1)
				v[it->node2] = (BlockId)i;
	}
	return v;
}

}