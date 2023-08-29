/****
DIAMOND protein aligner
Copyright (C) 2020 QIAGEN A/S (Aarhus, Denmark)
Code developed by Patrick Ettenhuber <patrick.ettenhuber@qiagen.com>

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

#pragma once
#include <string>
#include <vector>
#include "../basic/value.h"
#include "../data/sequence_file.h"
#include "../util/data_structures/flat_array.h"
#include "../basic/match.h"
#include "../dp/flags.h"
#include "../output/output_format.h"
#include "../util/algo/algo.h"

class ClusteringAlgorithm {
public:
	virtual void run() = 0;
	virtual std::string get_description() = 0;
	virtual ~ClusteringAlgorithm(){};
};

namespace Cluster {

struct CentroidSorted {};

void realign(const FlatArray<OId>& clusters, const std::vector<OId>& centroids, SequenceFile& db, std::function<void(const HspContext&)>& callback, HspValues hsp_values);
void realign(const std::vector<OId>& clustering, SequenceFile& db, std::function<void(const HspContext&)>& callback, HspValues hsp_values);
template<typename Int>
std::pair<FlatArray<Int>, std::vector<Int>> read(const std::string& file_name, const SequenceFile& db, CentroidSorted);
template<typename Int>
std::vector<Int> read(const std::string& file_name, const SequenceFile& db);
template<typename Int>
std::vector<Int> member2centroid_mapping(const FlatArray<Int>& clusters, const std::vector<Int>& centroids);
template<typename Int>
void output_mem(Util::Tsv::File& out, SequenceFile& db, const std::vector<Int>& mapping);
void output_mem(Util::Tsv::File& out, SequenceFile& db, Util::Tsv::File& oid_to_centroid_oid);
template<typename Int>
void output_mem(Util::Tsv::File& out, SequenceFile& db, std::vector<std::pair<Int, Int>>& mapping);
template<typename Int>
std::pair<FlatArray<Int>, std::vector<Int>> cluster_sorted(const std::vector<Int>& mapping);
template<typename Int>
std::pair<std::vector<Int>, std::vector<Int>> split(const std::vector<Int>& mapping);
std::vector<SuperBlockId> member_counts(const std::vector<SuperBlockId>& mapping);
Util::Tsv::File* open_out_tsv();
void init_thresholds();
std::vector<BlockId> len_sorted_clust(const FlatArray<Util::Algo::Edge<SuperBlockId>>& edges);

template<typename Int, typename Int2>
std::vector<Int2> convert_mapping(const std::vector<Int>& mapping, Int2) {
	std::vector<Int2> out;
	out.reserve(mapping.size());
	std::transform(mapping.begin(), mapping.end(), std::back_inserter(out), [](Int x) { return (Int2)x; });
	return out;
}

struct Mapback : public Consumer {
	Mapback(int64_t count) :
		centroid_id(count, -1),
		count(0)
	{}
	virtual void consume(const char* ptr, size_t n) {
		const char* end = ptr + n;
		OId query = -1;
		for(const char* p = ptr; p < end; p += sizeof(Output::Format::Edge::Data)) {
			const auto edge = *(Output::Format::Edge::Data*)p;
			assert(query == -1 || query == edge.query);
			query = edge.query;
			if (edge.qcovhsp >= config.member_cover)
				centroid_id[edge.query] = edge.target;
		}
		if (query >= 0 && centroid_id[query] == -1) {
			for (const char* p = ptr; p < end; p += sizeof(Output::Format::Edge::Data)) {
				const auto edge = *(Output::Format::Edge::Data*)p;
				if (edge.scovhsp >= config.member_cover) {
					covered_centroids.write((OId)query);
					covered_centroids.write((OId)edge.target);
					++count;
				}
			}
		}
	}
	std::vector<OId> unmapped() const {
		std::vector<OId> v;
		for (OId i = 0; i < (OId)centroid_id.size(); ++i)
			if (centroid_id[i] == -1)
				v.push_back(i);
		return v;
	}
	std::vector<std::pair<OId, OId>> targets_covered() {
		std::vector<std::pair<OId, OId>> v;
		v.resize(count);
		InputFile f(covered_centroids);
		f.read(v.data(), count);
		f.close_and_delete();
		return v;
	}
	std::vector<OId> centroid_id;
	TempFile covered_centroids;
	int64_t count;
};

template<typename It, typename It2>
int64_t update_clustering(It clustering, It2 mapping, It2 query_begin, It2 query_end, It2 db_begin) {
	const int64_t n = query_end - query_begin;
	int64_t k = 0;
	for (OId i = 0; i < n; ++i)
		if (mapping[i] >= 0 && clustering[query_begin[i]] != db_begin[mapping[i]]) {
			clustering[query_begin[i]] = (typename It::value_type)db_begin[mapping[i]];
			++k;
		}
	return k;
}

template<typename It>
std::vector<OId> cluster_members(It begin, It end, const FlatArray<OId>& clusters) {
	std::vector<OId> out;
	for (It i = begin; i != end; ++i) {
		for (auto it = clusters.cbegin(*i); it != clusters.cend(*i); ++it) {
			out.push_back(*it);
		}
	}
	return out;
}

}