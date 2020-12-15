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
#include <unordered_map>

using namespace std;

namespace Workflow { namespace Cluster {
	
string MultiStep::get_description() {
	return "A greedy stepwise vortex cover algorithm";
}

BitVector MultiStep::rep_bitset(const vector<int> &centroid, const BitVector *superset) {
	BitVector r(centroid.size());
	for (int c : centroid)
		if (!superset || superset->get(c))
			r.set(c);
	return r;
}

vector<int> MultiStep::cluster(DatabaseFile& db, const BitVector* filter) {
	statistics.reset();
	config.command = Config::blastp;
	//config.no_self_hits = true;
	//config.output_format = { "6", "qnum", "snum" };
	//config.output_format = { "6", "qnum", "snum", "qcovhsp", "scovhsp", "bitscore" };
	config.output_format = { "bin" };
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

	//message_stream << "Edges = " << nb.edges.size() << endl;

	unordered_map <uint32_t, NodEdgSet> components = find_connected_components(nb.smallest_index, nb.number_edges);
	message_stream << "Number of connected components: " << components.size() << endl;
	message_stream << "Average number of nodes per connected component: " << (double)nb.smallest_index.size() / components.size() << endl;

	uint32_t large = max_element(components.begin(), components.end(), [](const pair<uint32_t, NodEdgSet>& left, const pair<uint32_t, NodEdgSet>& right) {return left.second.nodes < right.second.nodes; })-> second.nodes;
	message_stream << "Largest connected component has " << large << " nodes." << endl;
	
	mapping_comp_set(components);

	uint32_t number_sets = max_element(components.begin(), components.end(), [](const pair<uint32_t, NodEdgSet>& left, const pair<uint32_t, NodEdgSet>& right) {return left.second.set < right.second.set; })->second.set;
	message_stream << "Number of sets: " << number_sets << endl;
	
 	return Util::Algo::greedy_vortex_cover(nb);

}

unordered_map<uint32_t, NodEdgSet> MultiStep::find_connected_components(vector<uint32_t>& sindex, const vector <size_t>& nedges){
	vector<uint32_t> seen_nodes;
	unordered_map<uint32_t, NodEdgSet> ne;

	uint32_t curr = 0;
	for (size_t k = 0; k < sindex.size(); k++) {
		curr = sindex[k];
		while (sindex[curr] != curr) {
			seen_nodes.push_back(sindex[curr]);
			curr = sindex[curr];
		}

		sindex[k] = curr;
		for (size_t i : seen_nodes) {
			sindex[i] = curr;
		}
		seen_nodes.clear();
	}


	for (size_t i = 0; i < sindex.size(); i++) {
		++ne[sindex[i]].nodes;
		ne[sindex[i]].edges += nedges[i];
	}
	return ne;
}

void MultiStep::mapping_comp_set(unordered_map<uint32_t, NodEdgSet>& comp) {
	vector <vector<uint32_t>> set;
	vector <size_t> size;
	
	bool TooBig;

	for (auto& it : comp) {
		TooBig = true;
		for (size_t j = 0; j < set.size(); j++) {
			if ((it.second.edges + size[j]) <= config.max_size_set) {
				set[j].push_back(it.first);
				size[j] += it.second.edges;
				comp[it.first].set = j;
				TooBig = false;
				break;
			}
		}
		if (TooBig || set.empty()) {
			set.push_back({ it.first });
			size.push_back({ it.second.edges });
			comp[it.first].set = set.size() - 1;
		}
	}
}


void MultiStep::steps(BitVector& current_reps, BitVector& previous_reps, vector <int>& current_centroids, vector <int>& previous_centroids, int count) {

	if (count == 0) {
		current_reps = rep_bitset(current_centroids);
	}
	else {
		current_reps = rep_bitset(current_centroids, &previous_reps);

		for (size_t i = 0; i < current_centroids.size(); ++i)
			if (!previous_reps.get(i))
				current_centroids[i] = current_centroids[previous_centroids[i]];
	}

	size_t n_rep_1 = previous_reps.one_count();
	size_t n_rep_2 = current_reps.one_count();

	if (n_rep_1 == 0) {
		n_rep_1 = current_centroids.size();
	}

	message_stream << "Clustering step " << count + 1 << " complete. #Input sequences: " << n_rep_1 << " #Clusters: " << n_rep_2 << endl;

	previous_centroids = move(current_centroids);
	previous_reps = move(current_reps);
}
		
void MultiStep::run() {
	if (config.database == "")
		throw runtime_error("Missing parameter: database file (--db/-d)");
	config.command = Config::makedb;
	unique_ptr<DatabaseFile> db(DatabaseFile::auto_create_from_fasta());
	const size_t seq_count = db->ref_header.sequences;

	BitVector current_reps;
	BitVector previous_reps;

	vector<int> current_centroids;
	vector<int> previous_centroids;
	
	for (size_t i = 0; i < config.cluster_steps.size(); i++) {
		
		if (config.sens_map.find(config.cluster_steps[i]) == config.sens_map.end()) {
			throw std::runtime_error("Invalid value for parameter --cluster-steps");
		}
			
		config.sensitivity = config.sens_map[config.cluster_steps[i]];
		current_centroids = cluster(*db, i==0 ? nullptr: &previous_reps);
		steps(current_reps, previous_reps, current_centroids, previous_centroids, i);
	}
		
	task_timer timer("Generating output");
	Sequence_set* rep_seqs;
	String_set<char, 0>* rep_ids;
	vector<unsigned> rep_database_id, rep_block_id(seq_count);
	db->rewind();
	db->load_seqs(&rep_database_id, (size_t)1e11, &rep_seqs, &rep_ids, true, &previous_reps);
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
		const unsigned r = rep_block_id[previous_centroids[i]];
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
