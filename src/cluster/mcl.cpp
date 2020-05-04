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
#include "mcl.h"

using namespace std;

namespace Workflow { namespace Cluster{
string MCL::get_key() {
	return "mcl";
}
string MCL::get_description() {
	return "Markov clustering according to doi:10.1137/040608635";
}

void MCL::get_exp(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r){
	if( r - (int) r == 0 ){
		// TODO: at some r it may be more beneficial to diagnoalize in and only take the exponents of the eigenvalues
		*out = *in;
		for(uint32_t i=1; i<r; i++){ 
			(*out) = (*out)*(*in);
		}
	}
	else{ 
		throw runtime_error(" Eigen does not provide an eigenvalue solver for sparse matrices");
	}
	*out = out->pruned();
}

void MCL::get_exp(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r){
	if( r - (int) r == 0 ){
		// TODO: at some r it may be more beneficial to diagnoalize in and only take the exponents of the eigenvalues
		*out = *in;
		for(uint32_t i=1; i<r; i++){ 
			(*out) *= (*in);
		}
	}
	else{ 
		Eigen::EigenSolver<Eigen::MatrixXf> solver(*in);
		Eigen::MatrixXcf d = solver.eigenvalues().asDiagonal();
		for(uint32_t idiag = 0; idiag<d.rows(); idiag++){
			d(idiag, idiag) = pow(d(idiag, idiag), r);
		}
		Eigen::MatrixXcf V = solver.eigenvectors();
		(*out) = (V * d * V.transpose()).real();
	}
}

void MCL::get_gamma(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r){
	vector<Eigen::Triplet<float>> data;
	for (uint32_t k=0; k<in->outerSize(); ++k){
		float colSum = 0.0f;
		for (Eigen::SparseMatrix<float>::InnerIterator it(*in, k); it; ++it){
			colSum += pow(it.value(),r);
		}
		for (Eigen::SparseMatrix<float>::InnerIterator it(*in, k); it; ++it){
			data.emplace_back(it.row(), it.col(), pow(it.value(),r) / colSum);
		}
	}
	out->setFromTriplets(data.begin(), data.end(), [] (const float&, const float &b) { return b; });
	*out = out->pruned(1.0, std::numeric_limits<float>::epsilon());
}

void MCL::get_gamma(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r){
	for (uint32_t icol=0; icol<in->cols(); ++icol){
		float colSum = 0.0f;
		for (uint32_t irow=0; irow<in->rows(); ++irow){
			colSum += pow((*in)(irow, icol) , r);
		}
		for (uint32_t irow=0; irow<in->rows(); ++irow){
			(*out)(irow, icol) =  pow((*in)(irow, icol) , r)/colSum;
		}
	}
}

vector<unordered_set<uint32_t>> MCL::get_list(Eigen::SparseMatrix<float>* m){
	LazyDisjointIntegralSet<uint32_t> disjointSet(m->cols());
	for (uint32_t k=0; k<m->outerSize(); ++k){
		for (Eigen::SparseMatrix<float>::InnerIterator it(*m, k); it; ++it){
			disjointSet.merge(it.row(), it.col());
		}
	}
	return disjointSet.getListOfSets();
}

vector<unordered_set<uint32_t>> MCL::get_list(Eigen::MatrixXf* m){
	LazyDisjointIntegralSet<uint32_t> disjointSet(m->cols());
	for (uint32_t icol=0; icol<m->cols(); ++icol){
		for (uint32_t irow=0; irow<m->rows(); ++irow){
			if( abs((*m)(irow, icol)) > 1e-7 ){
				disjointSet.merge(irow, icol);
			}
		}
	}
	return disjointSet.getListOfSets();
}

vector<unordered_set<uint32_t>> MCL::markov_process(Eigen::SparseMatrix<float>* m){
	uint32_t iteration = 0;
	double diff_norm = std::numeric_limits<double>::max();
	Eigen::SparseMatrix<float> msquared(m->rows(), m->cols());
	Eigen::SparseMatrix<float> m_update(m->rows(), m->cols());
	get_gamma(m, m, 1); // This is to get a matrix of random walks on the graph -> TODO: find out if something else is more suitable
	while( iteration < 100 && diff_norm > 1e-6*m->rows() ){
		get_exp(m, &msquared, 2);
		get_gamma(&msquared, &m_update, 2);
		*m -= m_update;
		diff_norm = m->norm();
		*m = m_update;
		iteration++;
	}
	return get_list(m);
}

vector<unordered_set<uint32_t>> MCL::markov_process(Eigen::MatrixXf* m){
	uint32_t iteration = 0;
	double diff_norm = std::numeric_limits<double>::max();
	Eigen::MatrixXf msquared(m->rows(), m->cols());
	Eigen::MatrixXf m_update(m->rows(), m->cols());
	get_gamma(m, m, 1); // This is to get a matrix of random walks on the graph -> TODO: find out if something else is more suitable
	while( iteration < 100 && diff_norm > 1e-6*m->rows() ){
		get_exp(m, &msquared, 2);
		get_gamma(&msquared, &m_update, 2);
		*m -= m_update;
		diff_norm = m->norm();
		*m = m_update;
		iteration++;
	}
	return get_list(m);
}

void MCL::run(){
	if (config.database == "")
		throw runtime_error("Missing parameter: database file (--db/-d)");
	config.command = Config::makedb;
	unique_ptr<DatabaseFile> db(DatabaseFile::auto_create_from_fasta());

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
	// 
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

	const size_t seq_count = db->ref_header.sequences;
	statistics.reset();
	config.command = Config::blastp;
	config.no_self_hits = false;
	config.output_format = { "6", "qnum", "snum", "qcovhsp", "scovhsp", "bitscore", "pident" };
	config.query_cover = 80;
	config.subject_cover = 80;
	config.algo = 0;
	config.index_mode = 0;
	config.freq_sd = 0;
	config.max_alignments = numeric_limits<uint64_t>::max();

	Workflow::Search::Options opt;
	opt.db = &(*db);
	opt.self = true;
	SparseMatrixStream<float> ms(db->ref_header.sequences);
	opt.consumer = &ms;
	opt.db_filter = nullptr;
	Workflow::Search::run(opt);

	task_timer timer;
	timer.go("Computing independent components");
	auto componentInfo = ms.getComponents();
	vector<vector<uint32_t>> indices = get<0>(componentInfo);
	vector<vector<Eigen::Triplet<float>>> components = get<1>(componentInfo);
	uint32_t nComponents = count_if(indices.begin(), indices.end(), [](vector<uint32_t> v){ return v.size() > 0;});
	uint32_t nComponentsLt1 = count_if(indices.begin(), indices.end(), [](vector<uint32_t> v){ return v.size() > 1;});
	timer.finish();
	cout << "Found " << nComponentsLt1 <<" ("<<nComponents<< " incl. singletons) disconnected components" << endl;

	timer.go("Clustering components");
	bool full = false;
	vector<uint32_t> clustering_result(db->ref_header.sequences, std::numeric_limits<uint32_t>::max());
	float max_sparsity = 0.0f;
	float min_sparsity = 1.0f;
	uint32_t cluster_id = 0;
	uint32_t nClustersEq1 = 0;
	uint32_t n_dense_calculations = 0;
	uint32_t n_sparse_calculations = 0;

	if( full ){
		// TODO: According to the SIAM publication this is not valid, just for debugging
		Eigen::SparseMatrix<float> m = ms.getMatrix();
		for(unordered_set<uint32_t> subset : markov_process(&m)){
			for(uint32_t el : subset){
				clustering_result[el] = cluster_id;
			}
			cluster_id++;
			if(subset.size() == 1) nClustersEq1 ++;
		}
	}
	else {
		// TODO: check if using dense matrices for non-sparse components improves performance
		for(uint32_t iComponent = 0; iComponent<indices.size(); iComponent++){
			vector<uint32_t> order = indices[iComponent];
			if(order.size() > 1){
				vector<Eigen::Triplet<float>> m = components[iComponent];
				assert(m.size() <= order.size() * order.size());
				float sparsity = 1.0-(1.0 * m.size()) / (order.size()*order.size());
				max_sparsity = max(max_sparsity, sparsity);
				min_sparsity = min(min_sparsity, sparsity);
				// map back to original ids
				vector<unordered_set<uint32_t>> list_of_sets;
				if(sparsity >= 0.8){ //TODO: a size limit should be given as well
					n_sparse_calculations++;
					Eigen::SparseMatrix<float> m_sparse(order.size(), order.size());
					m_sparse.setFromTriplets(m.begin(), m.end());
					list_of_sets = markov_process(&m_sparse);
				}
				else{
					n_dense_calculations++;
					Eigen::MatrixXf m_dense(order.size(), order.size());
					for(Eigen::Triplet<float> const & t : m){
						m_dense(t.row(), t.col()) = t.value();
					}
					list_of_sets = markov_process(&m_dense);
				}
				for(unordered_set<uint32_t> subset : list_of_sets){
					for(uint32_t el : subset){
						clustering_result[order[el]] = cluster_id;
					}
					cluster_id++;
					if(subset.size() == 1) nClustersEq1 ++;
				}
			}
			else if (order.size() == 1){
				clustering_result[order[0]] = cluster_id++;
				nClustersEq1++;
			}
		}
	}
	timer.finish();
	cout << "Found "<<cluster_id - nClustersEq1<< " ("<< cluster_id <<" incl. singletons) clusters with min sparsity "<<min_sparsity<< " and max. sparsity " << max_sparsity << " with " << n_dense_calculations << " dense, and " << n_sparse_calculations << " sparse calculations " << endl;

	timer.go("Cluster output");
	ostream *out = config.output_file.empty() ? &cout : new ofstream(config.output_file.c_str());
	vector<Letter> seq;
	string id;
	db->seek_direct();
	Hsp hsp;
	size_t n;
	out->precision(3);
	for (int i = 0; i < (int)db->ref_header.sequences; ++i) {
		db->read_seq(id, seq);
		(*out) << blast_id(id) << '\t' << clustering_result[i] << endl;		
		id.clear();
		seq.clear();
	}
	if (out != &cout) delete out;
	db->close();
	timer.finish();
}
}}