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

#include <sstream>
#include "cluster.h"
#include "util/io/file.h"
#include "util/string/tokenizer.h"
#include "basic/config.h"
#include "util/log_stream.h"
#include "util/system/system.h"

using std::pair;
using std::endl;
using std::vector;
using std::ostream;
using std::runtime_error;
using std::string;
using std::for_each;
using std::floor;
using std::unique_ptr;
using std::back_inserter;
using std::less;
using std::stringstream;

namespace Cluster {

const char* const HEADER_LINE = "centroid\tmember";

template<typename Int>
pair<FlatArray<Int>, vector<Int>> read(const string& file_name, const SequenceFile& db, CentroidSorted) {
	const int64_t lines = File::count_lines(file_name);
	File in(file_name, "rb", File::Flags::DETECT_COMPRESSION);
	string centroid, member;
	vector<pair<Int, Int>> pairs;
	pairs.reserve(lines);
	int64_t mappings = 0;
	const char* l;
	if (TabularFormat::header_format(::Config::cluster) == Header::SIMPLE) {
		l = in.getline();
		if (strcmp(l, HEADER_LINE) != 0)
			throw runtime_error("Clusters file is missing header line.");
	}
	while (l = in.getline(), !in.eof() || l[0] != '\0') {
		Util::String::Tokenizer<Util::String::CharDelimiter>(l, Util::String::CharDelimiter('\t')) >> centroid >> member;
		const Int centroid_oid = (Int)db.accession_to_oid(centroid).front();
		const Int member_oid = (Int)db.accession_to_oid(member).front();
		pairs.emplace_back(centroid_oid, member_oid);
		++mappings;
		if (mappings % 1000000 == 0)
			*log_stream << "#Entries: " << mappings << endl;
	}
	in.close();
	return make_flat_array(pairs.begin(), pairs.end(), config.threads_);
}

template pair<FlatArray<uint32_t>, vector<uint32_t>> read(const string&, const SequenceFile&, CentroidSorted);
template pair<FlatArray<uint64_t>, vector<uint64_t>> read(const string&, const SequenceFile&, CentroidSorted);

template<typename Int>
vector<Int> read(const string& file_name, const SequenceFile& db) {
	const int64_t lines = File::count_lines(file_name);
	File in(file_name, "rb", File::Flags::DETECT_COMPRESSION);
	string centroid, member;
	vector<Int> v(db.sequence_count().value());
	uint64_t mappings = 0;
	const char* l;
	if (TabularFormat::header_format(::Config::cluster) == Header::SIMPLE) {
		l = in.getline();
		if (strcmp(l, HEADER_LINE) != 0)
			throw runtime_error("Clustering input file is missing header line.");
	}
	while (l = in.getline(), !in.eof() || l[0] != '\0') {
		Util::String::Tokenizer<Util::String::CharDelimiter>(l, Util::String::CharDelimiter('\t')) >> centroid >> member;
		const auto centroid_oid = db.accession_to_oid(centroid);
		const auto member_oid = db.accession_to_oid(member);
		v[member_oid.front()] = (Int)centroid_oid.front();
		++mappings;
		if (mappings % 1000000 == 0)
			*log_stream << "#Entries: " << mappings << endl;
	}
	in.close();
	if (mappings != db.sequence_count().value())
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
	if (!config.approx_min_id.present() && config.min_id == 0.0)
		config.approx_min_id = config.command == ::Config::DEEPCLUST ? 0.0 : (config.command == ::Config::LINCLUST ? 90.0 : 50.0);
	//if (config.soft_masking.empty()) TODO
		//config.soft_masking = "tantan";
	// TODO
	return;
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

static int64_t seq_mem_use(Loc len, Loc id_len, int c, int min, int sketch_size) {
	assert(min > 1 || sketch_size > 0);
	int extend_stage = (min > 1 ? (len / (min / 2)) : sketch_size) * (15  // trace point buffer
		+ 16)        // SeedHitList
		+ 12         // SeedHitList
		+ 2 * len;   // SequenceSet
	extend_stage /= 2;
	return std::max(
		len + 8  // SequenceSet
		+ id_len + 8 // Seqtitle
		+ (min > 1 ? len / (min / 2) : sketch_size) * 9 / c // SeedArray
		+ 8 // super_block_id_to_oid
		+ 8 // BestCentroid / clustering
		+ 4 // unaligned 
		, extend_stage);
}

vector<string> cluster_steps(double id, bool linear) {
	if (!config.cluster_steps.empty()) {
		for (vector<string>::const_iterator it = config.cluster_steps.begin() + 1; it != config.cluster_steps.end(); ++it) {
			if (ends_with(*it, "_lin") && !ends_with(*(it - 1), "_lin"))
				throw runtime_error("Invalid cluster step sequence: linear steps must be at the start of the list.");
		}
		return config.cluster_steps;
	}
	vector<string> v = { "faster_lin" };
	if (id < 90)
		v.push_back("fast_lin");
	if (linear && id < 80) {
		if (id >= 70)
			v.push_back("linclust-70_lin");
		else if (id >= 60)
			v.push_back("linclust-60_lin");
		else if (id >= 40)
			v.push_back("linclust-40_lin");
		else
			v.push_back("linclust-20_lin");
	}
	else if (!linear) {
		if (id >= 80)
			v.push_back("fast");
		else
			v.push_back("default");
		if (id < 50)
			v.push_back("more-sensitive");
	}
	return v;
}

bool is_linclust(const vector<string>& steps) {
	for (const string& s : steps) {
		if (!ends_with(s, "_lin"))
			return false;
	}
	return true;
}

vector<string> default_round_approx_id(int steps) {
	return {};
	/*switch (steps) {
	case 1:
		return {};
	case 2:
		return { "27.0" };
	default:
		return { "27.0", "0.0" };
	}*/
}

vector<string> default_round_cov(int steps) {
	return {};
	/*switch (steps) {
	case 1:
		return {};
	case 2:
		return { "85.0" };
	default:
		return { "87.0", "85.0" };
	}*/
}

static int round_ccd(const string& depth) {
	stringstream ss(depth);
	int i;
	ss >> i;
	if (ss.fail() || !ss.eof())
		throw runtime_error("Invalid number format for --connected-component-depth");
	return i;
}

int round_ccd(std::vector<std::string> param, int round, int round_count, bool linear) {
	if (param.size() > 1 && param.size() != (size_t)round_count)
		throw runtime_error("Parameter count for --connected-component-depth has to be 1 or the number of cascaded clustering rounds.");
	if (param.size() == 0)
		return 0;
	if (param.size() > 1)
		return round_ccd(param[round]);
	return (round == round_count - 1) ^ linear ? 1 : round_ccd(param[0]);
}

}