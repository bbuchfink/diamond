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

vector<int> MultiStep::cluster(SequenceFile& db, const BitVector* filter) {
	statistics.reset();
	config.command = Config::blastp;
	//config.no_self_hits = true;
	//config.output_format = { "6", "qnum", "snum" };
	//config.output_format = { "6", "qnum", "snum", "qcovhsp", "scovhsp", "bitscore" };
	config.output_format = { "bin" };
	config.query_cover = 80;
	config.subject_cover = 80;
	config.algo = Config::Algo::DOUBLE_INDEXED;
	config.freq_sd = 0;
	config.max_alignments = numeric_limits<size_t>::max();

	Workflow::Search::Options opt;
	opt.db = &db;
	opt.self = true;
	Neighbors nb(db.sequence_count());

	opt.consumer = &nb;
	opt.db_filter = filter;

	Workflow::Search::run(opt);
	
	/*
	auto lo = nb.dSet.getListOfSets();
	for (auto& set : lo) {
		cout << "set :";
		for (auto item : set) {
			cout << " " << item;
		}
		cout << endl;
	}
	*/

	vector<unordered_set<uint32_t>> connected = nb.dSet.getListOfSets();
	vector<uint32_t> EdgSet(nb.number_edges.size());
	unordered_map <uint32_t, NodEdgSet> components = find_connected_components(connected, EdgSet, nb.number_edges);
	message_stream << "Number of connected components: " << components.size() << endl;
	message_stream << "Average number of nodes per connected component: " << (double)nb.number_edges.size() / components.size() << endl;

	uint32_t large = max_element(components.begin(), components.end(), [](const pair<uint32_t, NodEdgSet>& left, const pair<uint32_t, NodEdgSet>& right) {return left.second.nodes < right.second.nodes; })-> second.nodes;
	message_stream << "Largest connected component has " << large << " nodes." << endl;

	vector<TempFile*> tmp_sets = mapping_comp_set(components);

	uint32_t number_sets = max_element(components.begin(), components.end(), [](const pair<uint32_t, NodEdgSet>& left, const pair<uint32_t, NodEdgSet>& right) {return left.second.set < right.second.set; })-> second.set;
	message_stream << "Number of sets: " << number_sets + 1 << endl;


	if (config.external) {
		save_edges_external(nb.tempfiles, tmp_sets, components, EdgSet);
		return cluster_sets(db.sequence_count(), tmp_sets);
	}

	else {
		return Util::Algo::greedy_vertex_cover(nb);
	}

}

void MultiStep::save_edges_external(vector<TempFile*>& all_edges, vector<TempFile*>& sorted_edges, const unordered_map<uint32_t, NodEdgSet>& comp, const vector<uint32_t>& s_index){

	size_t count = 0;
	uint32_t query;
	uint32_t subject;
	uint32_t result_set;
	while (count < all_edges.size()) {
		InputFile tmp_edges(*all_edges[count]);
		while (true) {
			try {
				tmp_edges.read(query);
				tmp_edges.read(subject);		
				result_set = comp.at(s_index[query]).set;
				sorted_edges[result_set]->write(query);
				sorted_edges[result_set]->write(subject);
			}
			catch (EndOfStream&) {
				tmp_edges.close_and_delete();
				break;
			}
		}
		count++;
	}
}

vector<int> MultiStep::cluster_sets(const size_t nb_size, vector<TempFile*> &sorted_edges){
	vector<int> cluster(nb_size);
	vector<vector<int>> tmp_neighbors(nb_size);
	vector<int>curr;
	uint32_t query;
	uint32_t subject;

	iota(cluster.begin(), cluster.end(), 0);

	for (size_t i = 0; i < sorted_edges.size(); i++) {
		InputFile tmp(*sorted_edges[i]);
		delete(sorted_edges[i]);
		while (true) {
			try {
				tmp.read(query);
				tmp.read(subject);
				tmp_neighbors[query].push_back(subject);
			}
			catch (EndOfStream&) {
				tmp.close_and_delete();
				break;
			}
		}
		curr = Util::Algo::greedy_vertex_cover(tmp_neighbors);

		for (int i = 0; i < (int)curr.size(); i++) {
				if (curr[i] != i) {
					cluster[i] = curr[i];
				}
		}

		tmp_neighbors.clear();
		tmp_neighbors.resize(nb_size);
	}

	return cluster;
}

unordered_map<uint32_t, NodEdgSet> MultiStep::find_connected_components(const vector<unordered_set<uint32_t>> &connected, vector<uint32_t> &EdgSet, const vector <size_t>& nedges){
	
	unordered_map<uint32_t, NodEdgSet> ne;
	for (size_t i = 0; i < connected.size(); i++) {
		for(uint32_t j : connected[i]){
			EdgSet[j] = i;
			++ne[i].nodes; 
			ne[i].edges += nedges[j];
		}
	}

	return ne;
}

vector<TempFile*> MultiStep::mapping_comp_set(unordered_map<uint32_t, NodEdgSet>& comp) {
	vector <vector<uint32_t>> set;
	vector <size_t> size;
	vector<TempFile*> temp_set;
	
	
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
			if (config.external) {
				temp_set.push_back(new TempFile());
			}
		}
	}
	return temp_set;
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
	unique_ptr<SequenceFile> db(SequenceFile::auto_create());
	const size_t seq_count = db->sequence_count();

	BitVector current_reps;
	BitVector previous_reps;

	vector<int> current_centroids;
	vector<int> previous_centroids;
	
	for (size_t i = 0; i < config.cluster_steps.size(); i++) {
		
		config.sensitivity = from_string<Sensitivity>(config.cluster_steps[i]);
		current_centroids = cluster(*db, i==0 ? nullptr: &previous_reps);
		steps(current_reps, previous_reps, current_centroids, previous_centroids, i);
	}
		
	task_timer timer("Generating output");
	SequenceSet* rep_seqs;
	String_set<char, 0>* rep_ids;
	vector<unsigned> rep_database_id, rep_block_id(seq_count);
	db->set_seqinfo_ptr(0);
	db->load_seqs(&rep_database_id, (size_t)1e11, &rep_seqs, &rep_ids, true, &previous_reps);
	for (size_t i = 0; i < rep_database_id.size(); ++i)
		rep_block_id[rep_database_id[i]] = (unsigned)i;

	ostream* out = config.output_file.empty() ? &cout : new ofstream(config.output_file.c_str());
	vector<Letter> seq;
	string id;
	db->init_seq_access();
	Hsp hsp;
	out->precision(3);

	for (int i = 0; i < (int)db->sequence_count(); ++i) {
		db->read_seq(seq, id);
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
