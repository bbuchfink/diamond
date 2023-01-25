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

#include "../basic/config.h"
#include "../util/log_stream.h"
#include "cluster.h"
#include "../basic/statistics.h"
#include "../run/workflow.h"
#include "../data/fasta/fasta_file.h"
#include "cascaded/cascaded.h"

using std::endl;
using std::shared_ptr;
using std::for_each;
using std::vector;
using std::tie;
using std::make_shared;
using std::unique_ptr;

namespace Cluster {

void reassign() {
	config.database.require();
	config.clustering.require();
	message_stream << "Coverage cutoff: " << config.member_cover << '%' << endl;

	task_timer timer("Opening the database");
	shared_ptr<SequenceFile> db(SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NEED_LETTER_COUNT | SequenceFile::Flags::ACC_TO_OID_MAPPING | SequenceFile::Flags::OID_TO_ACC_MAPPING, SequenceFile::Metadata()));
	config.db_size = db->letters();
	timer.finish();
	message_stream << "#Database sequences: " << db->sequence_count() << ", #Letters: " << db->letters() << endl;
	unique_ptr<Util::Tsv::File> out(open_out_tsv());

	timer.go("Reading the input file");
	vector<OId> clustering = read<OId>(config.clustering, *db);

	timer.go("Finding centroids");
	vector<OId> centroids, members;
	tie(centroids, members) = split(clustering);

	timer.go("Creating member database");
	shared_ptr<FastaFile> member_db(db->sub_db(members.cbegin(), members.cend()));
	member_db->set_seqinfo_ptr(0);

	timer.go("Creating centroid database");
	shared_ptr<FastaFile> centroid_db(db->sub_db(centroids.cbegin(), centroids.cend()));
	centroid_db->set_seqinfo_ptr(0);
	timer.finish();

	statistics.reset();
	init_thresholds();
	config.command = Config::blastp;
	config.max_target_seqs_ = 1;
	config.output_format = { "edge" };
	config.self = false;
	config.query_cover = config.member_cover;
	config.sensitivity = from_string<Sensitivity>(cluster_steps(config.approx_min_id).back());
	shared_ptr<Mapback> mapback = make_shared<Mapback>(members.size());
	Search::run(centroid_db, member_db, mapback);

	timer.go("Updating clustering");
	const int64_t n = update_clustering(clustering.begin(), mapback->centroid_id.cbegin(), members.cbegin(), members.cend(), centroids.cbegin());
	timer.finish();
	
	message_stream << "Reassigned members: " << n << '/' << members.size() << endl;

	timer.go("Generating output");
	if (flag_any(db->format_flags(), SequenceFile::FormatFlags::TITLES_LAZY))
		db->init_random_access(0, 0, false);
	output_mem(*out, *db, clustering);

	timer.go("Closing the database");
	db.reset();
}

}