/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include "multi_step_cluster.h"

using namespace std;

namespace Workflow { namespace Cluster{ 
string MultiStep::get_key() {
	return "multi-step";
}
string MultiStep::get_description(){
	return "A greedy stepwise vortex cover algorithm";
}

vector<bool> MultiStep::rep_bitset(const vector<int> &centroid, const vector<bool> *superset) {
	vector<bool> r(centroid.size());
	for (int c : centroid)
		if(!superset || (*superset)[c])
			r[c] = true;
	return r;
}

vector<int> MultiStep::cluster(DatabaseFile &db, const vector<bool> *filter) {
	statistics.reset();
	config.command = Config::blastp;
	config.no_self_hits = true;
	//config.output_format = { "6", "qnum", "snum" };
	config.output_format = { "6", "qnum", "snum", "qcovhsp", "scovhsp", "bitscore" };
	config.query_cover = 80;
	config.subject_cover = 80;
	config.algo = 0;
	config.index_mode = 0;
	config.freq_sd = 0;
	config.max_alignments = numeric_limits<size_t>::max();

	Workflow::Search::Options opt;
	opt.db = &db;
	opt.self = true;
	Neighbors nb(db.ref_header.sequences);
	opt.consumer = &nb;
	opt.db_filter = filter;

	Workflow::Search::run(opt);

	return Util::Algo::greedy_vortex_cover(nb);
}

void MultiStep::run() {
	if (config.database == "")
		throw runtime_error("Missing parameter: database file (--db/-d)");
	config.command = Config::makedb;
	unique_ptr<DatabaseFile> db(DatabaseFile::auto_create_from_fasta());
	const size_t seq_count = db->ref_header.sequences;

	config.min_id = 70;
	vector<int> centroid1(cluster(*db, nullptr));
	vector<bool> rep1(rep_bitset(centroid1));
	size_t n_rep1 = count_if(rep1.begin(), rep1.end(), [](bool x) { return x; });
	message_stream << "Clustering step 1 complete. #Input sequences: " << centroid1.size() << " #Clusters: " << n_rep1 << endl;

	config.mode_more_sensitive = true;
	config.min_id = 0;
	vector<int> centroid2(cluster(*db, &rep1));
	vector<bool> rep2(rep_bitset(centroid2, &rep1));
	for (size_t i = 0; i < centroid2.size(); ++i)
		if (!rep1[i])
			centroid2[i] = centroid2[centroid1[i]];
	message_stream << "Clustering step 2 complete. #Input sequences: " << n_rep1 << " #Clusters: " << count_if(rep2.begin(), rep2.end(), [](bool x) { return x; }) << endl;

	task_timer timer("Generating output");
	Sequence_set *rep_seqs;
	String_set<char, 0> *rep_ids;
	vector<unsigned> rep_database_id, rep_block_id(seq_count);
	db->rewind();
	db->load_seqs(&rep_database_id, (size_t)1e11, &rep_seqs, &rep_ids, true, &rep2);
	for (size_t i = 0; i < rep_database_id.size(); ++i)
		rep_block_id[rep_database_id[i]] = (unsigned)i;

	ostream *out = config.output_file.empty() ? &cout : new ofstream(config.output_file.c_str());
	vector<Letter> seq;
	string id;
	db->seek_direct();
	Hsp hsp;
	size_t n;
	out->precision(3);

	for (int i = 0; i < (int)db->ref_header.sequences; ++i) {
		db->read_seq(id, seq);
		const unsigned r = rep_block_id[centroid2[i]];
		(*out) << blast_id(id) << '\t'
			<< blast_id((*rep_ids)[r]) << '\t';		
		if ((int)i == centroid2[i])
			(*out) << "100\t100\t100\t0" << endl;
		else {
			Masking::get().bit_to_hard_mask(seq.data(), seq.size(), n);
			smith_waterman(sequence(seq), (*rep_seqs)[r], hsp);
			(*out) << hsp.id_percent() << '\t'
				<< hsp.query_cover_percent((unsigned)seq.size()) << '\t'
				<< hsp.subject_cover_percent((unsigned)(*rep_seqs)[r].length()) << '\t'
				<< score_matrix.bitscore(hsp.score) << endl;
		}
	}
	if (out != &cout) delete out;
	delete rep_seqs;
	delete rep_ids;
	db->close();
}
}}
