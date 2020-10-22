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

namespace Workflow { namespace Cluster {
	
string MultiStep::get_key() {
	return "multi-step";
}
string MultiStep::get_description() {
	return "A greedy stepwise vortex cover algorithm";
}

vector<bool> MultiStep::rep_bitset(const vector<int>& centroid, const vector<bool>* superset) {
	vector<bool> r(centroid.size());
	for (int c : centroid)
		if (!superset || (*superset)[c])
			r[c] = true;
	return r;
}

vector<int> MultiStep::cluster(DatabaseFile& db, const vector<bool>* filter) {
	statistics.reset();
	config.command = Config::blastp;
	//config.no_self_hits = true;
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

	message_stream << "Edges = " << nb.edges.size() << endl;

	return Util::Algo::greedy_vortex_cover(nb);
}


void MultiStep::steps(vector <vector<bool>> &reps, vector <vector <int>> &centroids, int count) {
	count = count + 1;

	if (reps.empty()) {
		vector <bool> rep_c(rep_bitset(centroids.front()));
		reps.push_back(rep_c);
	}

	else {
		vector<bool> rep_c(rep_bitset(centroids.back(), &reps.back()));
		reps.push_back(rep_c);

		for (size_t i = 0; i < centroids.back().size(); ++i)
			if (!reps.front()[i])
				centroids.back()[i] = centroids.back()[centroids.front()[i]];
	}

	size_t n_rep_1 = count_if(reps.front().begin(), reps.front().end(), [](bool x) { return x; });
	size_t n_rep_2 = count_if(reps.back().begin(), reps.back().end(), [](bool x) { return x; });

	if (n_rep_1 == n_rep_2) {
		message_stream << "Clustering step " << count << " complete. #Input sequences: " << centroids.front().size() << " #Clusters: " << n_rep_1 << endl;
	}

	else {
		message_stream << "Clustering step " << count << " complete. #Input sequences: " << n_rep_1 << " #Clusters: " << n_rep_2 << endl;
		
		centroids.front().swap(centroids.back());
		centroids.pop_back();

		reps.front().swap(reps.back());
		reps.pop_back();
	}
}
		
void MultiStep::run() {
	if (config.database == "")
		throw runtime_error("Missing parameter: database file (--db/-d)");
	config.command = Config::makedb;
	unique_ptr<DatabaseFile> db(DatabaseFile::auto_create_from_fasta());
	const size_t seq_count = db->ref_header.sequences;

	vector<vector<int>> centroids;
	vector<vector<bool>> reps;

	
	for (int i = 0; i <= config.cluster_steps.size(); i++) {
		int c = i - 1;

		if (i == 0) {
			vector<int> centroid_current(cluster(*db, nullptr));
			centroids.push_back(centroid_current);
		}

		else {
			if (config.cluster_steps[c] != "fast" && config.cluster_steps[c] != "sensitive" && config.cluster_steps[c] != "mid-sensitive" &&
				config.cluster_steps[c] != "more-sensitive" && config.cluster_steps[c] != "very-sensitive" &&
				config.cluster_steps[c] != "ultra-sensitive") {
				throw std::runtime_error("Invalid value for parameter --cluster-steps");
			}

			else {
				config.sensitivity = config.sens_map[config.cluster_steps[c]];
				vector<int> centroid_current(cluster(*db, &reps.back()));
				centroids.push_back(centroid_current);
			}
		}
		steps(reps, centroids, i);
	}
	

	task_timer timer("Generating output");
	Sequence_set* rep_seqs;
	String_set<char, 0>* rep_ids;
	vector<unsigned> rep_database_id, rep_block_id(seq_count);
	db->rewind();
	db->load_seqs(rep_database_id, (size_t)1e11, &rep_seqs, &rep_ids, true, &reps.back());
	for (size_t i = 0; i < rep_database_id.size(); ++i)
		rep_block_id[rep_database_id[i]] = (unsigned)i;

	ostream* out = config.output_file.empty() ? &cout : new ofstream(config.output_file.c_str());
	vector<Letter> seq;
	string id;
	db->seek_direct();
	Hsp hsp;
	size_t n;
	out->precision(3);

	for (int i = 0; i < (int)db->ref_header.sequences; ++i) {
		db->read_seq(id, seq);
		const unsigned r = rep_block_id[centroids.back()[i]];
		(*out) << blast_id(id) << '\t'
			<< blast_id((*rep_ids)[r]) << '\n';
		/*if ((int)i == centroid2[i])
			(*out) << "100\t100\t100\t0" << endl;
		else {
			Masking::get().bit_to_hard_mask(seq.data(), seq.size(), n);
			smith_waterman(sequence(seq), (*rep_seqs)[r], hsp);
			(*out) << hsp.id_percent() << '\t'
				<< hsp.query_cover_percent((unsigned)seq.size()) << '\t'
				<< hsp.subject_cover_percent((unsigned)(*rep_seqs)[r].length()) << '\t'
				<< score_matrix.bitscore(hsp.score) << endl;
		} */
	}
	if (out != &cout) delete out;
	delete rep_seqs;
	delete rep_ids;
	db->close();

}
}}