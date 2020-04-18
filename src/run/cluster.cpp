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

#include <Eigen/Sparse>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <stdio.h>
#include <memory>
#include <fstream>
#include <limits>
#include "../util/system/system.h"
#include "../util/util.h"
#include "../basic/config.h"
#include "../data/reference.h"
#include "workflow.h"
#include "disjoint_set.h"
#include "../util/io/consumer.h"
#include "../util/algo/algo.h"
#include "../basic/statistics.h"
#include "../util/log_stream.h"
#include "../dp/dp.h"
#include "../basic/masking.h"

using namespace std;

namespace Workflow { namespace Cluster {

struct Neighbors : public vector<vector<int>>, public Consumer {
	Neighbors(size_t n):
		vector<vector<int>>(n)
	{}
	virtual void consume(const char *ptr, size_t n) override {
		int query, subject, count;
		float qcov, scov, bitscore;
		const char *end = ptr + n;
		while (ptr < end) {
			//if (sscanf(ptr, "%i\t%i\n%n", &query, &subject, &count) != 2)
			if (sscanf(ptr, "%i\t%i\t%f\t%f\t%f\n%n", &query, &subject, &qcov, &scov, &bitscore, &count) != 5)
				throw runtime_error("Cluster format error.");
			ptr += count;
			//cout << query << '\t' << subject << '\t' << qcov << '\t' << scov << '\t' << endl;
			(*this)[query].push_back(subject);
			edges.push_back({ query, subject, (int)bitscore });
		}
	}
	vector<Util::Algo::Edge> edges;
};

class NeighborStream : public Consumer {
	LazyDisjointSet<size_t>* disjointSet;
	virtual void consume(const char *ptr, size_t n) override {
		size_t query, subject, count;
		float qcov, scov, bitscore;
		const char *end = ptr + n;
		while (ptr < end) {
			//if (sscanf(ptr, "%i\t%i\n%n", &query, &subject, &count) != 2)
			if (sscanf(ptr, "%lu\t%lu\t%f\t%f\t%f\n%ln", &query, &subject, &qcov, &scov, &bitscore, &count) != 5)
				throw runtime_error("Cluster format error.");
			disjointSet->merge(query, subject);
			ptr += count;
		}
	}
public:
	NeighborStream(size_t n){
		disjointSet = new LazyDisjointIntegralSet<size_t>(n); 
	}
	~NeighborStream(){
		delete disjointSet;
	}
	vector<unordered_set<size_t>> getListOfSets(){
		return disjointSet->getListOfSets();
	}
};

template <typename T>
class SparseMatrixStream : public Consumer {
	size_t n;
	vector<Eigen::Triplet<T>> data;
	LazyDisjointSet<size_t>* disjointSet;
	virtual void consume(const char *ptr, size_t n) override {
		size_t query, subject, count;
		float qcov, scov, bitscore, id;
		const char *end = ptr + n;
		while (ptr < end) {
			//if (sscanf(ptr, "%i\t%i\n%n", &query, &subject, &count) != 2)
			if (sscanf(ptr, "%lu\t%lu\t%f\t%f\t%f\t%f\n%ln", &query, &subject, &qcov, &scov, &bitscore, &id, &count) != 6) 
				throw runtime_error("Cluster format error.");
			data.emplace_back(query, subject, (qcov/100.0f) * (scov/100.0f) * (id/100.0f));
			disjointSet->merge(query, subject);
			ptr += count;
		}
	}
public:
	SparseMatrixStream(size_t n){
		this->n = n;
		disjointSet = new LazyDisjointIntegralSet<size_t>(n); 
		// Note: this makes sure the self hits are always included, the value needs to be aligned with the measure used in consume!
		for(uint32_t idiag=0; idiag<n; idiag++){
			data.emplace_back(idiag, idiag, 1.0);
		}
	}
	~SparseMatrixStream(){
		delete disjointSet;
	}
	pair<vector<vector<uint32_t>>, vector<Eigen::SparseMatrix<T>>> getComponents(){
		vector<unordered_set<size_t>> sets = disjointSet->getListOfSets();
		vector<vector<Eigen::Triplet<T>>> split(sets.size());
		unordered_map<size_t, uint32_t> indexToSetId;
		for(uint32_t iset = 0; iset < sets.size(); iset++){
			split.push_back(vector<Eigen::Triplet<T>>());
			for(size_t index : sets[iset]){
				indexToSetId.emplace(index, iset);
			}
		}
		for(Eigen::Triplet<T> d : data){
			assert(indexToSetId[d.row()] == indexToSetId[d.row()]);
			split[indexToSetId[d.row()]].emplace_back(d.row(), d.col(), d.value());
		}
		vector<vector<uint32_t>> indices(sets.size());
		vector<Eigen::SparseMatrix<T>> components(sets.size());
		for(uint32_t iset = 0; iset < sets.size(); iset++){
			Eigen::SparseMatrix<T> component(sets[iset].size(), sets[iset].size());
			vector<uint32_t> order(sets[iset].begin(), sets[iset].end());
			map<uint32_t, uint32_t> index_map;
			uint32_t iel = 0;
			for(uint32_t el: order){
				index_map.emplace(el, iel++);
			}
			vector<Eigen::Triplet<T>> remapped(split[iset].size());
			for(Eigen::Triplet<T> t : split[iset]){
				remapped.emplace_back(index_map[t.row()], index_map[t.col()], t.value());
			}
			component.setFromTriplets(remapped.begin(), remapped.end());
			components.emplace_back(move(component));
			indices.emplace_back(move(order));
		}
		return std::make_pair(indices, components);
	}
	Eigen::SparseMatrix<T> getMatrix(){
		Eigen::SparseMatrix<T> m(n,n);
		m.setFromTriplets(data.begin(), data.end());
		return m;
	}
};

vector<bool> rep_bitset(const vector<int> &centroid, const vector<bool> *superset = nullptr) {
	vector<bool> r(centroid.size());
	for (int c : centroid)
		if(!superset || (*superset)[c])
			r[c] = true;
	return r;
}

vector<int> cluster(DatabaseFile &db, const vector<bool> *filter) {
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
	config.max_alignments = numeric_limits<uint64_t>::max();

	Workflow::Search::Options opt;
	opt.db = &db;
	opt.self = true;
	Neighbors nb(db.ref_header.sequences);
	opt.consumer = &nb;
	opt.db_filter = filter;

	Workflow::Search::run(opt);

	return Util::Algo::greedy_vortex_cover(nb);
}

void run_two_step_clustering(DatabaseFile& db){
	const size_t seq_count = db.ref_header.sequences;

	config.min_id = 70;
	vector<int> centroid1(cluster(db, nullptr));
	vector<bool> rep1(rep_bitset(centroid1));
	size_t n_rep1 = count_if(rep1.begin(), rep1.end(), [](bool x) { return x; });
	message_stream << "Clustering step 1 complete. #Input sequences: " << centroid1.size() << " #Clusters: " << n_rep1 << endl;

	config.mode_more_sensitive = true;
	config.min_id = 0;
	vector<int> centroid2(cluster(db, &rep1));
	vector<bool> rep2(rep_bitset(centroid2, &rep1));
	for (size_t i = 0; i < centroid2.size(); ++i)
		if (!rep1[i])
			centroid2[i] = centroid2[centroid1[i]];
	message_stream << "Clustering step 2 complete. #Input sequences: " << n_rep1 << " #Clusters: " << count_if(rep2.begin(), rep2.end(), [](bool x) { return x; }) << endl;

	task_timer timer("Generating output");
	Sequence_set *rep_seqs;
	String_set<char, 0> *rep_ids;
	vector<unsigned> rep_database_id, rep_block_id(seq_count);
	db.rewind();
	db.load_seqs(rep_database_id, (size_t)1e11, &rep_seqs, &rep_ids, true, &rep2);
	for (size_t i = 0; i < rep_database_id.size(); ++i)
		rep_block_id[rep_database_id[i]] = (unsigned)i;

	ostream *out = config.output_file.empty() ? &cout : new ofstream(config.output_file.c_str());
	vector<Letter> seq;
	string id;
	db.seek_direct();
	Hsp hsp;
	size_t n;
	out->precision(3);

	for (int i = 0; i < (int)db.ref_header.sequences; ++i) {
		db.read_seq(id, seq);
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
}


void run_transitive_closure_clustering(DatabaseFile& db){
	// testing
	// LazyDisjointIntegralSet<uint32_t> test(5);
	// test.merge(1,3);
	// test.merge(3,4);
	// test.merge(0,2);
	// auto los = test.getListOfSets();
	// for(auto& set: los){
	// 	for(auto item:set){
	// 		cout<< " " << item;
	// 	}
	// 	cout<<endl;
	// }
	// unordered_set<string> s ( {"red","green","blue","rain", "sun", "snow"} );
	// LazyDisjointTypeSet<string> ts(&s);

	// ts.merge("red", "green");
	// ts.merge("snow", "rain");
	// ts.merge("sun", "rain");
	// auto lo = ts.getListOfSets();
	// for(auto& set: lo){
	// 	cout<< "set :";
	// 	for(auto item:set){
	// 		cout<< " " << item;
	// 	}
	// 	cout<<endl;
	// }

	const size_t seq_count = db.ref_header.sequences;
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
	config.max_alignments = numeric_limits<uint64_t>::max();

	Workflow::Search::Options opt;
	opt.db = &db;
	opt.self = true;
	NeighborStream ns(db.ref_header.sequences);
	opt.consumer = &ns;
	opt.db_filter = nullptr;

	Workflow::Search::run(opt);

	auto clusterResult = ns.getListOfSets();
	cout << "Found "<<clusterResult.size()<<" clusters"<<endl;
	// for(auto& set: clusterResult){
	// 	cout<< "set :";
	// 	for(auto item : set){
	// 		cout<< " " << item;
	// 	}
	// 	cout<<endl;
	// }
}

Eigen::SparseMatrix<float> get_gamma(Eigen::SparseMatrix<float>* m, float r){
	vector<Eigen::Triplet<float>> data;
	for (uint32_t k=0; k<m->outerSize(); ++k){
		float colSum = 0.0f;
		for (Eigen::SparseMatrix<float>::InnerIterator it(*m, k); it; ++it){
			colSum += pow(it.value(),r);
		}
		for (Eigen::SparseMatrix<float>::InnerIterator it(*m, k); it; ++it){
			data.emplace_back(it.row(), it.col(), pow(it.value(),r) / colSum);
		}
	}
	Eigen::SparseMatrix<float> gamma(m->rows(), m->cols());
	gamma.setFromTriplets(data.begin(), data.end());
	return gamma.pruned(1.0, std::numeric_limits<float>::epsilon());
}

vector<unordered_set<uint32_t>> get_list(Eigen::SparseMatrix<float>* m){
	LazyDisjointIntegralSet<uint32_t> disjointSet(m->cols());
	for (uint32_t k=0; k<m->outerSize(); ++k){
		for (Eigen::SparseMatrix<float>::InnerIterator it(*m, k); it; ++it){
			disjointSet.merge(it.row(), it.col());
		}
	}
	return disjointSet.getListOfSets();
}

vector<unordered_set<uint32_t>> markov_process(Eigen::SparseMatrix<float>* m){
	uint32_t iteration = 0;
	double diff_norm = std::numeric_limits<double>::max();
	*m  = get_gamma(m,1); // This is to get a matrix of random walks on the graph -> TODO: find out if something else is more suitable
	while( iteration < 100 && diff_norm > 1e-6*m->rows() ){
		Eigen::SparseMatrix<float> msquared = ((*m) * (*m)).pruned(1.0, std::numeric_limits<float>::epsilon());
		Eigen::SparseMatrix<float> m_update = get_gamma(&msquared, 2);
		diff_norm = (*m - m_update).norm();
		*m = m_update.pruned(1.0, std::numeric_limits<float>::epsilon());
		iteration++;
	}
	return get_list(m);
}

void run_markov_clustering(DatabaseFile& db){
	// testing
	// vector<Eigen::Triplet<float>> data;
	// data.emplace_back(0, 0, 0.2f);
	// data.emplace_back(0, 1, 0.25f);
	// data.emplace_back(0, 5, 0.333f);
	// data.emplace_back(0, 6, 0.25f);
	// data.emplace_back(0, 9, 0.25f);
	// data.emplace_back(1, 0, 0.20f);
	// data.emplace_back(1, 1, 0.25f);
	// data.emplace_back(1, 2, 0.25f);
	// data.emplace_back(1, 4, 0.20f);
	// data.emplace_back(2, 1, 0.25f);
	// data.emplace_back(2, 2, 0.25f);
	// data.emplace_back(2, 3, 0.20f);
	// data.emplace_back(2, 4, 0.20f);
	// data.emplace_back(3, 2, 0.25f);
	// data.emplace_back(3, 3, 0.20f);
	// data.emplace_back(3, 7, 0.20f);
	// data.emplace_back(3, 8, 0.20f);
	// data.emplace_back(3, 10, 0.20f);
	// data.emplace_back(4, 1, 0.25f);
	// data.emplace_back(4, 2, 0.25f);
	// data.emplace_back(4, 4, 0.20f);
	// data.emplace_back(4, 6, 0.25f);
	// data.emplace_back(4, 7, 0.20f);
	// data.emplace_back(5, 0, 0.20f);
	// data.emplace_back(5, 5, 0.333f);
	// data.emplace_back(5, 9, 0.25f);
	// data.emplace_back(6, 0, 0.20f);
	// data.emplace_back(6, 4, 0.20f);
	// data.emplace_back(6, 6, 0.25f);
	// data.emplace_back(6, 9, 0.25f);
	// data.emplace_back(7, 3, 0.20f);
	// data.emplace_back(7, 4, 0.20f);
	// data.emplace_back(7, 7, 0.20f);
	// data.emplace_back(7, 8, 0.20f);
	// data.emplace_back(7, 10, 0.20f);
	// data.emplace_back(8, 3, 0.20f);
	// data.emplace_back(8, 7, 0.20f);
	// data.emplace_back(8, 8, 0.20f);
	// data.emplace_back(8, 10, 0.20f);
	// data.emplace_back(8, 11, 0.333f);
	// data.emplace_back(9, 0, 0.20f);
	// data.emplace_back(9, 5, 0.333f);
	// data.emplace_back(9, 6, 0.25f);
	// data.emplace_back(9, 9, 0.25f);
	// data.emplace_back(10, 3, 0.20f);
	// data.emplace_back(10, 7, 0.20f);
	// data.emplace_back(10, 8, 0.20f);
	// data.emplace_back(10, 10, 0.20f);
	// data.emplace_back(10, 11, 0.333f);
	// data.emplace_back(11, 8, 0.20f);
	// data.emplace_back(11, 10, 0.20f);
	// data.emplace_back(11, 11, 0.333f);
	// Eigen::SparseMatrix<float> tmp(12,12);
	// tmp.setFromTriplets(data.begin(), data.end());
	// auto cr = markov_process(&tmp);
	// cout << "Found "<<cr.size()<<" clusters"<<endl;
	// for(auto& set: cr){
	// 	cout<< "set :";
	// 	for(auto item : set){
	// 		cout<< " " << item;
	// 	}
	// 	cout<<endl;
	// }
	// return;

	const size_t seq_count = db.ref_header.sequences;
	statistics.reset();
	config.command = Config::blastp;
	config.no_self_hits = true;
	config.output_format = { "6", "qnum", "snum", "qcovhsp", "scovhsp", "bitscore", "pident" };
	config.query_cover = 80;
	config.subject_cover = 80;
	config.algo = 0;
	config.index_mode = 0;
	config.freq_sd = 0;
	config.max_alignments = numeric_limits<uint64_t>::max();

	Workflow::Search::Options opt;
	opt.db = &db;
	opt.self = true;
	SparseMatrixStream<float> ms(db.ref_header.sequences);
	opt.consumer = &ms;
	opt.db_filter = nullptr;
	Workflow::Search::run(opt);

	auto componentInfo = ms.getComponents();
	vector<vector<uint32_t>> indices = get<0>(componentInfo);
	vector<Eigen::SparseMatrix<float>> components = get<1>(componentInfo);
	uint32_t nComponents = count_if(indices.begin(), indices.end(), [](vector<uint32_t> v){ return v.size() > 0;});
	uint32_t nComponentsLt1 = count_if(indices.begin(), indices.end(), [](vector<uint32_t> v){ return v.size() > 1;});

	cout << "DIAMOND done" << endl;
	cout << "************" << endl;
	cout << "Found " << nComponentsLt1 <<" ("<<nComponents<< " incl. singletons) disconnected components" << endl;

	vector<unordered_set<uint32_t>> clustering_result;
	float max_sparsity = 0.0f;
	float min_sparsity = 1.0f;
	bool full = false;
	if( full ){
		// TODO: According to the SIAM publication this is not valid, just for debugging
		Eigen::SparseMatrix<float> m = ms.getMatrix();
		for(unordered_set<uint32_t> subset : markov_process(&m)){
			clustering_result.emplace_back(subset);
		}
	}
	else {
		// TODO: check if using dense matrices for non-sparse components improves performance
		for(uint32_t iComponent = 0; iComponent<indices.size(); iComponent++){
			vector<uint32_t> order = indices[iComponent];
			if(order.size() > 1){
				Eigen::SparseMatrix<float> m = components[iComponent];
				float sparsity = 1.0-(1.0 * m.nonZeros()) / (m.rows()*m.cols());
				max_sparsity = max(max_sparsity, sparsity);
				min_sparsity = min(min_sparsity, sparsity);
				// map back to original ids
				for(unordered_set<uint32_t> subset : markov_process(&m)){
					unordered_set<uint32_t> s;
					for(uint32_t el : subset){
						s.emplace(order[el]);
					}
					clustering_result.emplace_back(s);
				}
			}
			else if (order.size() == 1){
				clustering_result.emplace_back(order.begin(), order.end());
			}
		}
	}
	uint32_t nClusters = clustering_result.size();
	uint32_t nClustersLt1 = count_if(clustering_result.begin(), clustering_result.end(), [](unordered_set<uint32_t> v){ return v.size() > 1;});
	cout << "Found "<<nClustersLt1<< " ("<< nClusters<<" incl. singletons) clusters with min sparsity "<<min_sparsity<< " and max. sparsity " << max_sparsity << endl;
	// for(auto& set: clustering_result){
	// 	cout<< "set :";
	// 	for(auto item : set){
	// 		cout<< " " << item;
	// 	}
	// 	cout<<endl;
	// }
}

void run() {
	if (config.database == "")
		throw runtime_error("Missing parameter: database file (--db/-d)");
	config.command = Config::makedb;
	unique_ptr<DatabaseFile> db(DatabaseFile::auto_create_from_fasta());
	switch(config.cluster_algo){
		case Config::multi_step:
			run_two_step_clustering(*db);
			break;
		case Config::transitive_closure:
			run_transitive_closure_clustering(*db);
			break;
		case Config::markov_clustering:
			run_markov_clustering(*db);
			break;
	}
	db->close();
}
}}
