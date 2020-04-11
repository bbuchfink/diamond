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
#include "sparse_matrix.h"
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
	SparseSimpleMatrix<T>* m = nullptr;
	LazyDisjointSet<size_t>* disjointSet;
	virtual void consume(const char *ptr, size_t n) override {
		size_t query, subject, count;
		float qcov, scov, bitscore;
		const char *end = ptr + n;
		while (ptr < end) {
			//if (sscanf(ptr, "%i\t%i\n%n", &query, &subject, &count) != 2)
			if (sscanf(ptr, "%lu\t%lu\t%f\t%f\t%f\n%ln", &query, &subject, &qcov, &scov, &bitscore, &count) != 5) 
				throw runtime_error("Cluster format error.");
			m->set_elm(query, subject, bitscore);
			disjointSet->merge(query, subject);
			ptr += count;
		}
	}
public:
	SparseMatrixStream(size_t n, bool symmetric){
		disjointSet = new LazyDisjointIntegralSet<size_t>(n); 
		if( symmetric ){
			m = new SparseSymmetricSimpleMatrix<T>(n);
		}
		else{
			m = new SparseSimpleMatrixImpl<T>(n, n);
		}
	}
	~SparseMatrixStream(){
		if(m) delete m;
		delete disjointSet;
	}
	SparseSimpleMatrix<T>* getMatrix(){
		cout << " getting matrix " << m->get_nrows() << " " << m->get_ncols() << " " << m->get_n_nonzero_elements() << " sparsity "<< (1.0 * m->get_n_nonzero_elements()) / (m->get_nrows()* m->get_ncols()) << endl;
		SparseSimpleMatrix<T>* h = m;
		m = nullptr;
		return h;
	}
	vector<unordered_set<size_t>> getListOfSets(){
		return disjointSet->getListOfSets();
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

SparseSimpleMatrix<float>* get_gamma(SparseSimpleMatrix<float>* m, uint32_t r){
	SparseSimpleMatrix<float>* gamma = new SparseSimpleMatrixImpl<float>(m->get_nrows(), m->get_ncols());
	float* column_sum = new float[m->get_ncols()];
	for(uint32_t icol = 0; icol< m->get_ncols(); icol++){
		column_sum[icol] = 0.0f;
	}
	for(auto it = m->begin(); it != m->end(); it++){
		uint32_t i,j;
		m->get_indices(it->first, &i, &j);
		float el = pow(it->second, r);
		gamma->set_elm(i,j,el);
		column_sum[j] += el;
	}
	for(auto it = gamma->begin(); it != gamma->end(); it++){
		uint32_t i,j;
		gamma->get_indices(it->first, &i, &j);
		gamma->set_elm(i,j, it->second/column_sum[j]);
	}
	return gamma;
}
SimpleMatrix<float>* get_gamma(SimpleMatrix<float>* m, uint32_t r){
	SimpleMatrix<float>* gamma = new SparseSimpleMatrixImpl<float>(m->get_nrows(), m->get_ncols());
	float* column_sum = new float[m->get_ncols()];
	for(uint32_t icol = 0; icol< m->get_ncols(); icol++){
		column_sum[icol] = 0.0f;
	}
	for(uint32_t irow = 0; irow < m->get_nrows(); irow++){
		for(uint32_t icol = 0; icol < m->get_ncols(); icol++){
			column_sum[icol]+=pow(m->get_elm(irow, icol), r);
		}
	}
	for(uint32_t irow = 0; irow < m->get_nrows(); irow++){
		for(uint32_t icol = 0; icol < m->get_ncols(); icol++){
			gamma->set_elm(irow, icol, pow(m->get_elm(irow, icol),r) / column_sum[icol]);
		}
	}
	return gamma;
}

SimpleMatrix<float>* get_exp(SimpleMatrix<float>* m, uint32_t r){
	// TODO: make this more efficient
	SimpleMatrix<float>* tmp1 = m->multiply(m);
	for(uint32_t i = 2; i<r; i++){
		SimpleMatrix<float>* tmp2 = tmp1->multiply(tmp1);
		delete tmp1;
		tmp1 = tmp2;
	}
	return tmp1;
}
vector<unordered_set<uint32_t>> get_list(SimpleMatrix<float>* m){
	LazyDisjointIntegralSet<uint32_t> disjointSet(m->get_ncols());
	for(uint32_t irow = 0; irow < m->get_nrows(); irow++){
		for(uint32_t icol = 0; icol < m->get_ncols(); icol++){
			if(abs(m->get_elm(irow, icol)) > 1e-13){
				disjointSet.merge(irow, icol);
			}
		}
	}
	return disjointSet.getListOfSets();
}
vector<unordered_set<uint32_t>> get_list(SparseSimpleMatrix<float>* m){
	LazyDisjointIntegralSet<uint32_t> disjointSet(m->get_ncols());
	for(auto it = m->begin(); it != m->end(); it++){
		uint32_t irow, icol;
		m->get_indices(it->first, &irow, &icol);
		disjointSet.merge(irow, icol);
	}
	return disjointSet.getListOfSets();
}


vector<unordered_set<uint32_t>> markov_process(SimpleMatrix<float>* m){
	uint32_t iteration = 0;
	double diff_norm = std::numeric_limits<double>::max();
	while( iteration < 100 && diff_norm > 1e-3 ){
		SimpleMatrix<float>* msquared = get_exp(m, 2);
		SimpleMatrix<float>* m_update = get_gamma(msquared, 2);
		delete msquared;
		SimpleMatrix<float>* diff_m = m_update->minus(m);
		diff_norm = diff_m->norm(2.0,2.0);
		delete m;
		delete diff_m;
		m = m_update;
		m_update = nullptr;
		iteration++;
	}
	cout<< "Markov Process converged after " << iteration <<" iterations " << diff_norm << endl;
	return get_list(m);
}

void run_markov_clustering(DatabaseFile& db){
	// testing
	// SparseSimpleMatrix<float>* tmp = new SparseSimpleMatrixImpl<float>(12, 12);
	// tmp->set_elm(0, 0, 0.2f);
	// tmp->set_elm(0, 1, 0.25f);
	// tmp->set_elm(0, 5, 0.333f);
	// tmp->set_elm(0, 6, 0.25f);
	// tmp->set_elm(0, 9, 0.25f);
	// tmp->set_elm(1, 0, 0.20f);
	// tmp->set_elm(1, 1, 0.25f);
	// tmp->set_elm(1, 2, 0.25f);
	// tmp->set_elm(1, 4, 0.20f);
	// tmp->set_elm(2, 1, 0.25f);
	// tmp->set_elm(2, 2, 0.25f);
	// tmp->set_elm(2, 3, 0.20f);
	// tmp->set_elm(2, 4, 0.20f);
	// tmp->set_elm(3, 2, 0.25f);
	// tmp->set_elm(3, 3, 0.20f);
	// tmp->set_elm(3, 7, 0.20f);
	// tmp->set_elm(3, 8, 0.20f);
	// tmp->set_elm(3, 10, 0.20f);
	// tmp->set_elm(4, 1, 0.25f);
	// tmp->set_elm(4, 2, 0.25f);
	// tmp->set_elm(4, 4, 0.20f);
	// tmp->set_elm(4, 6, 0.25f);
	// tmp->set_elm(4, 7, 0.20f);
	// tmp->set_elm(5, 0, 0.20f);
	// tmp->set_elm(5, 5, 0.333f);
	// tmp->set_elm(5, 9, 0.25f);
	// tmp->set_elm(6, 0, 0.20f);
	// tmp->set_elm(6, 4, 0.20f);
	// tmp->set_elm(6, 6, 0.25f);
	// tmp->set_elm(6, 9, 0.25f);
	// tmp->set_elm(7, 3, 0.20f);
	// tmp->set_elm(7, 4, 0.20f);
	// tmp->set_elm(7, 7, 0.20f);
	// tmp->set_elm(7, 8, 0.20f);
	// tmp->set_elm(7, 10, 0.20f);
	// tmp->set_elm(8, 3, 0.20f);
	// tmp->set_elm(8, 7, 0.20f);
	// tmp->set_elm(8, 8, 0.20f);
	// tmp->set_elm(8, 10, 0.20f);
	// tmp->set_elm(8, 11, 0.333f);
	// tmp->set_elm(9, 0, 0.20f);
	// tmp->set_elm(9, 5, 0.333f);
	// tmp->set_elm(9, 6, 0.25f);
	// tmp->set_elm(9, 9, 0.25f);
	// tmp->set_elm(10, 3, 0.20f);
	// tmp->set_elm(10, 7, 0.20f);
	// tmp->set_elm(10, 8, 0.20f);
	// tmp->set_elm(10, 10, 0.20f);
	// tmp->set_elm(10, 11, 0.333f);
	// tmp->set_elm(11, 8, 0.20f);
	// tmp->set_elm(11, 10, 0.20f);
	// tmp->set_elm(11, 11, 0.333f);
	// auto cr = markov_process(tmp);
	// cout << "Found "<<cr.size()<<" clusters"<<endl;
	// for(auto& set: cr){
	// 	cout<< "set :";
	// 	for(auto item : set){
	// 		cout<< " " << item;
	// 	}
	// 	cout<<endl;
	// }
	// delete tmp;
	// return;


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
	SparseMatrixStream<float> ms(db.ref_header.sequences, false);
	opt.consumer = &ms;
	opt.db_filter = nullptr;
	Workflow::Search::run(opt);

	SparseSimpleMatrix<float>* m = ms.getMatrix();

	auto clusterResult = markov_process(m);
	cout << "Found "<<clusterResult.size()<<" clusters"<<endl;
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
