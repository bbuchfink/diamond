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
	std::chrono::high_resolution_clock::time_point t = std::chrono::high_resolution_clock::now();
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
	sparse_exp_time += (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count() / 1000.0;
}

void MCL::get_exp(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r){
	// TODO: at some r it may be more beneficial to diagnoalize in and only take the exponents of the eigenvalues
	std::chrono::high_resolution_clock::time_point t = std::chrono::high_resolution_clock::now();
	if(r - (int) r == 0){
		*out = *in;
		for(uint32_t i=1; i<r; i++){ 
			(*out) *= (*in);
		}
		dense_int_exp_time += (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count() / 1000.0;
	}
	else{ 
		// TODO: check whether the matrix is self-adjoint and use SelfAdjointEigenSolver instead
		// TODO: try and get out->noalias() = in->pow(r); to work (unsupported! http://eigen.tuxfamily.org/dox/unsupported/group__MatrixFunctions__Module.html).
		Eigen::EigenSolver<Eigen::MatrixXf> solver(*in);
		Eigen::MatrixXcf d = solver.eigenvalues().asDiagonal();
		for(uint32_t idiag = 0; idiag<d.rows(); idiag++){
			d(idiag, idiag) = pow(d(idiag, idiag), r);
		}
		Eigen::MatrixXcf V = solver.eigenvectors();
		double thr = 0.5 * abs(V.determinant());
		if( thr > std::numeric_limits<float>::epsilon() ){
			out->noalias() = (V * d.real() * V.inverse()).real();
		}
		dense_gen_exp_time += (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count() / 1000.0;
	}
}

void MCL::get_gamma(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r){
	std::chrono::high_resolution_clock::time_point t = std::chrono::high_resolution_clock::now();
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
	sparse_gamma_time += (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count() / 1000.0;
}

void MCL::get_gamma(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r){
	std::chrono::high_resolution_clock::time_point t = std::chrono::high_resolution_clock::now();
	// Note that Eigen matrices are column-major, so this is the most efficient way
	for (uint32_t icol=0; icol<in->cols(); ++icol){
		float colSum = 0.0f;
		for (uint32_t irow=0; irow<in->rows(); ++irow){
			colSum += pow(in->coeffRef(irow, icol) , r);
		}
		for (uint32_t irow=0; irow<in->rows(); ++irow){
			out->coeffRef(irow, icol) =  pow(in->coeffRef(irow, icol) , r) / colSum;
		}
	}
	dense_gamma_time += (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count() / 1000.0;
}

vector<unordered_set<uint32_t>> MCL::get_list(Eigen::SparseMatrix<float>* m){
	std::chrono::high_resolution_clock::time_point t = std::chrono::high_resolution_clock::now();
	LazyDisjointIntegralSet<uint32_t> disjointSet(m->cols());
	for (uint32_t k=0; k<m->outerSize(); ++k){
		for (Eigen::SparseMatrix<float>::InnerIterator it(*m, k); it; ++it){
			disjointSet.merge(it.row(), it.col());
		}
	}
	sparse_list_time += (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count() / 1000.0;
	return disjointSet.getListOfSets();
}

vector<unordered_set<uint32_t>> MCL::get_list(Eigen::MatrixXf* m){
	std::chrono::high_resolution_clock::time_point t = std::chrono::high_resolution_clock::now();
	LazyDisjointIntegralSet<uint32_t> disjointSet(m->cols());
	for (uint32_t icol=0; icol<m->cols(); ++icol){
		for (uint32_t irow=0; irow<m->rows(); ++irow){
			if( abs((*m)(irow, icol)) > std::numeric_limits<float>::epsilon()){
				disjointSet.merge(irow, icol);
			}
		}
	}
	dense_list_time += (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count() / 1000.0;
	return disjointSet.getListOfSets();
}

// void printMatrix(Eigen::SparseMatrix<float>* m){
// 	printf("[");
// 	for (uint32_t irow=0; irow<m->rows(); ++irow){
// 		printf("[");
// 		for (uint32_t icol=0; icol<m->cols(); ++icol){
// 			printf("%12.7f",m->coeffRef(irow, icol));
// 		}
// 		if(irow<m->rows() - 1){
// 			printf("];\n");
// 		}
// 		else{
// 			printf("]");
// 		}
// 	}
// 	printf("]\n\n");
// }
// void printMatrix(Eigen::MatrixXf* m){
// 	printf("[");
// 	for (uint32_t irow=0; irow<m->rows(); ++irow){
// 		printf("[");
// 		for (uint32_t icol=0; icol<m->cols(); ++icol){
// 			printf("%12.7f",m->coeffRef(irow, icol));
// 		}
// 		if(irow<m->rows() - 1){
// 			printf("];\n");
// 		}
// 		else{
// 			printf("]");
// 		}
// 	}
// 	printf("]\n\n");
// }

vector<unordered_set<uint32_t>> MCL::markov_process(Eigen::SparseMatrix<float>* m, float inflation, float expansion){
	for (uint32_t idiag=0; idiag<m->cols(); ++idiag){
		assert(abs(m->coeffRef(idiag, idiag)) > std::numeric_limits<float>::epsilon());
	}
	uint32_t iteration = 0;
	float diff_norm = std::numeric_limits<float>::max();
	Eigen::SparseMatrix<float> msquared(m->rows(), m->cols());
	Eigen::SparseMatrix<float> m_update(m->rows(), m->cols());
	get_gamma(m, m, 1); // This is to get a matrix of random walks on the graph -> TODO: find out if something else is more suitable
	while( iteration < 100 && diff_norm > 1e-6*m->rows() ){
		get_exp(m, &msquared, expansion);
		get_gamma(&msquared, &m_update, inflation);
		*m -= m_update;
		diff_norm = m->norm();
		*m = m_update;
		iteration++;
	}
	return get_list(m);
}

vector<unordered_set<uint32_t>> MCL::markov_process(Eigen::MatrixXf* m, float inflation, float expansion){
	for (uint32_t idiag=0; idiag<m->cols(); ++idiag){
		assert(abs(m->coeffRef(idiag, idiag)) > std::numeric_limits<float>::epsilon());
	}
	uint32_t iteration = 0;
	float diff_norm = std::numeric_limits<float>::max();
	Eigen::MatrixXf msquared(m->rows(), m->cols());
	Eigen::MatrixXf m_update(m->rows(), m->cols());
	get_gamma(m, m, 1); // This is to get a matrix of random walks on the graph -> TODO: find out if something else is more suitable
	while( iteration < 100 && diff_norm > std::numeric_limits<float>::epsilon()*m->rows() ){
		get_exp(m, &msquared, expansion);
		get_gamma(&msquared, &m_update, inflation);
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
	// vector<unordered_set<uint32_t>> cr;
	// float inf = 2.0;
	// float exp = 2.0;
	// if(false){
	// 	Eigen::SparseMatrix<float> tmp(12,12);
	// 	tmp.setFromTriplets(data.begin(), data.end());
	// 	printMatrix(&tmp);
	// 	cr = markov_process(&tmp, inf, exp);
	// 	printMatrix(&tmp);
	// }
	// else{
	// 	Eigen::MatrixXf tmp = Eigen::MatrixXf::Zero(12,12);
	// 	for(Eigen::Triplet<float> const & t : data){
	// 		tmp(t.row(), t.col()) = t.value();
	// 	}
	// 	printMatrix(&tmp);
	// 	cr = markov_process(&tmp, inf, exp);
	// 	printMatrix(&tmp);
	// }
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
	string format = config.cluster_similarity;
	if(format.empty()){
		format = "qcovhsp/100*scovhsp/100*pident/100";
	}
	config.output_format = {"clus", format};
	config.algo = 0;
	config.index_mode = 0;
	config.freq_sd = 0;

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
	std::vector<uint32_t> sort_order(components.size());
	std::iota(sort_order.begin(), sort_order.end(), 0);
	std::sort(sort_order.begin(), sort_order.end(), [&](uint32_t i, uint32_t j){return indices[i].size() > indices[j].size();});

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
	float inflation = (float) config.cluster_mcl_inflation;
	float expansion = (float) config.cluster_mcl_expansion;

	if( full ){
		// TODO: According to the SIAM publication this is not valid, just for debugging
		Eigen::SparseMatrix<float> m = ms.getMatrix();
		for(unordered_set<uint32_t> subset : markov_process(&m, inflation, expansion)){
			for(uint32_t el : subset){
				clustering_result[el] = cluster_id;
			}
			cluster_id++;
			if(subset.size() == 1) nClustersEq1 ++;
		}
	}
	else {
		// TODO: check when/if using dense matrices for non-sparse components improves performance
		// vector<std::thread> threads;
		// for (size_t i = 0; i < config.threads_; ++i)
		// 	threads.emplace_back(seed_join_worker, query_idx, ref_idx, &seedp, &range, query_seed_hits, ref_seed_hits);
		// for (auto &t : threads)
		// 	t.join();
		for(uint32_t iComponent : sort_order){
			vector<uint32_t> order = indices[iComponent];
			if(order.size() > 1){
				vector<Eigen::Triplet<float>> m = components[iComponent];
				assert(m.size() <= order.size() * order.size());
				float sparsity = 1.0-(1.0 * m.size()) / (order.size()*order.size());
				max_sparsity = max(max_sparsity, sparsity);
				min_sparsity = min(min_sparsity, sparsity);
				// map back to original ids
				vector<unordered_set<uint32_t>> list_of_sets;

				//TODO: a size limit for the dense matrix should control this as well
				std::chrono::high_resolution_clock::time_point t = std::chrono::high_resolution_clock::now();
				if(sparsity >= config.cluster_mcl_sparsity_switch && expansion - (int) expansion == 0){ 
					n_sparse_calculations++;
					Eigen::SparseMatrix<float> m_sparse(order.size(), order.size());
					m_sparse.setFromTriplets(m.begin(), m.end());
					sparse_create_time += (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count() / 1000.0;
					list_of_sets = markov_process(&m_sparse, inflation, expansion);
				}
				else{
					n_dense_calculations++;
					Eigen::MatrixXf m_dense = Eigen::MatrixXf::Zero(order.size(), order.size());
					for(Eigen::Triplet<float> const & t : m){
						m_dense(t.row(), t.col()) = t.value();
					}
					dense_create_time += (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count() / 1000.0;
					list_of_sets = markov_process(&m_dense, inflation, expansion);
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

	cout << "Time used for matrix creation: " << sparse_create_time + dense_create_time <<" (sparse: " << sparse_create_time << ", dense: " << dense_create_time <<")" << endl;
	cout << "Time used for exp: " << sparse_exp_time + dense_int_exp_time + dense_gen_exp_time << " (sparse: " << sparse_exp_time << ", dense int: " << dense_int_exp_time << ", dense gen: " << dense_gen_exp_time <<")" << endl; 
	cout << "Time used for gamma: " << sparse_gamma_time + dense_gamma_time <<" (sparse: " << sparse_gamma_time << ", dense: " << dense_gamma_time <<")" << endl;
	cout << "Time used for listing: " << sparse_list_time + dense_list_time <<" (sparse: " << sparse_list_time << ", dense: " << dense_list_time <<")" << endl;
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
