/****
DIAMOND protein aligner
Copyright (C) 2022 Max Planck Society for the Advancement of Science e.V.

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

#include "cascaded.h"
#include "../../basic/config.h"
#include "../../data/fasta/fasta_file.h"
#include "../../run/workflow.h"
#include "../../output/output_format.h"
#include "../../basic/statistics.h"
#include "../search/search.h"

using std::unique_ptr;
using std::endl;
using std::vector;
using std::function;
using std::make_shared;
using std::shared_ptr;
using std::for_each;
using std::to_string;
using std::string;
using std::tie;
using std::pair;

namespace Cluster {

static vector<OId> recluster(shared_ptr<SequenceFile>& db, const vector<OId>& clustering, int iteration) {
	task_timer timer(("*** Initializing recluster iteration " + to_string(iteration + 1)).c_str());

	FlatArray<OId> clusters;
	vector<OId> centroids;
	tie(clusters, centroids) = cluster_sorted(clustering);

	BitVector centroid_aligned(db->sequence_count());
	for_each(centroids.begin(), centroids.end(), [&centroid_aligned](OId c) { centroid_aligned.set(c); });
	function<void(const HspContext&)> callback([&centroid_aligned](const HspContext& h) {
		if (h.scovhsp() >= config.member_cover && (h.approx_id() >= config.approx_min_id || h.id_percent() >= config.approx_min_id))
			centroid_aligned.set(h.subject_oid);
	});
	timer.finish();

	HspValues hsp_values = HspValues::TARGET_COORDS;
	if (config.approx_min_id > 0)
		hsp_values |= HspValues::QUERY_COORDS | HspValues::IDENT | HspValues::LENGTH;
	realign(clusters, centroids, *db, callback, hsp_values);

	timer.go("Creating database of unaligned sequences");
	const vector<OId> unal_members = centroid_aligned.negative_list();
	if (unal_members.empty())
		return clustering;
	shared_ptr<FastaFile> unaligned(db->sub_db(unal_members.cbegin(), unal_members.cend()));
	unaligned->set_seqinfo_ptr(0);
	timer.finish();
	message_stream << "#Sequences that failed to align against assigned centroid: " << unal_members.size() << endl;

	timer.go("Creating centroid database");
	shared_ptr<FastaFile> centroid_db(db->sub_db(centroids.cbegin(), centroids.cend()));
	centroid_db->set_seqinfo_ptr(0);
	timer.finish();

	statistics.reset();
	config.command = Config::blastp;
	config.max_target_seqs_ = 1;
	config.iterate = vector<string>();
	config.output_format = { "edge" };
	config.self = false;
	config.query_cover = config.no_recluster_bd ? config.member_cover : 0;
	config.subject_cover = 0;
	config.query_or_target_cover = config.no_recluster_bd ? 0 : config.member_cover;
	config.sensitivity = from_string<Sensitivity>(cluster_steps(config.approx_min_id).back());
	//tie(config.chunk_size, config.lowmem_) = block_size(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)), Search::iterated_sens.at(config.sensitivity).front(), false);
	config.lowmem_ = 1;
	config.chunk_size = 4.0;
	shared_ptr<Mapback> mapback = make_shared<Mapback>(unal_members.size());
	Search::run(centroid_db, unaligned, mapback);

	timer.go("Updating clustering");
	vector<OId> out = clustering;
	update_clustering(out.begin(), mapback->centroid_id.cbegin(), unal_members.cbegin(), unal_members.cend(), centroids.cbegin());

	timer.go("Deallocating memory");	
	centroid_db.reset();
	timer.finish();

	vector<OId> unmapped_members = mapback->unmapped();
	message_stream << "#Sequences that failed to align against any centroid: " << unmapped_members.size() << endl;
	if (unmapped_members.empty())
		return out;

	shared_ptr<SequenceFile> unmapped;
	if (config.no_recluster_bd) {
		timer.go("Creating database of unmapped sequences");
		unmapped.reset(unaligned->sub_db(unmapped_members.cbegin(), unmapped_members.cend()));
	}
	else {
		timer.go("Loading covered centroids list");
		vector<pair<OId, OId>> covered_centroids = mapback->targets_covered();
		timer.finish();
		message_stream << "#Centroid sequences tentatively covered: " << covered_centroids.size() << endl;
		timer.go("Filtering list");
		ips4o::parallel::sort(covered_centroids.begin(), covered_centroids.end());
		vector<OId> centroid_list;
		auto it = covered_centroids.begin();
		auto it2 = unmapped_members.begin();
		while (it < covered_centroids.end() && it2 < unmapped_members.end()) {
			if (it->first < *it2)
				++it;
			else if (*it2 < it->first)
				++it2;
			else {
				centroid_list.push_back(it->second);
				++it;
			}
		}
		ips4o::parallel::sort(centroid_list.begin(), centroid_list.end());
		auto end = std::unique(centroid_list.begin(), centroid_list.end());
		const int64_t n = end - centroid_list.begin();
		timer.finish();
		message_stream << "#Centroid sequences covered: " << n << endl;
		timer.go("Making sequence list for reclustering");
		for (int64_t i = 0; i < unmapped_members.size(); ++i)
			unmapped_members[i] = unal_members[unmapped_members[i]];
		const vector<OId> members = cluster_members(centroid_list.begin(), end, clusters);
		unmapped_members.reserve(unmapped_members.size() + members.size());
		unmapped_members.insert(unmapped_members.end(), members.begin(), members.end());
		//unmapped_members.reserve(unmapped_members.size() + n);
		//for (auto i = centroid_list.begin(); i < end; ++i)
			//unmapped_members.push_back(centroids[*i]);
		ips4o::parallel::sort(unmapped_members.begin(), unmapped_members.end());
		timer.finish();
		message_stream << "#Sequences from covered clustes: " << members.size() << endl;
		timer.go("Creating database for reclustering");
		unmapped.reset(db->sub_db(unmapped_members.cbegin(), unmapped_members.cend()));
	}
	mapback.reset();
	timer.finish();

	timer.go("Deallocating memory");
	clusters.clear();
	clusters.shrink_to_fit();
	unaligned.reset();
	timer.finish();

	const vector<OId> reclust = recluster(unmapped, convert_mapping(cascaded(unmapped), OId()), iteration + 1);

	timer.go("Deallocating memory");
	unmapped.reset();	

	timer.go("Merging clusterings");
	for (SuperBlockId i = 0; i < (SuperBlockId)reclust.size(); ++i)
		out[unal_members[unmapped_members[i]]] = unal_members[unmapped_members[reclust[i]]];

	return out;
}

void recluster() {
	config.database.require();
	config.clustering.require();
	init_thresholds();
	//config.strict_gvc = true;
	message_stream << "Coverage cutoff: " << config.member_cover << '%' << endl;

	task_timer timer("Opening the database");
	shared_ptr<SequenceFile> db(SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NEED_LETTER_COUNT | SequenceFile::Flags::ACC_TO_OID_MAPPING | SequenceFile::Flags::OID_TO_ACC_MAPPING, SequenceFile::Metadata()));
	config.db_size = db->letters();
	timer.finish();
	unique_ptr<Util::Tsv::File> out(open_out_tsv());
	message_stream << "#Database sequences: " << db->sequence_count() << ", #Letters: " << db->letters() << endl;

	timer.go("Reading the input file");
	const vector<OId> clustering = read<OId>(config.clustering, *db);
	timer.finish();

	const vector<OId> member2centroid = recluster(db, clustering, 0);

	timer.go("Generating output");
	if (flag_any(db->format_flags(), SequenceFile::FormatFlags::TITLES_LAZY))
		db->init_random_access(0, 0, false);
	output_mem(*out, *db, member2centroid);

	timer.go("Closing the database");
	db.reset();
}

}