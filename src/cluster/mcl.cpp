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

#include <numeric>
#include "mcl.h"

#define MASK_INVERSE        0xC000000000000000
#define MASK_NORMAL_NODE    0x4000000000000000
#define MASK_ATTRACTOR_NODE 0x8000000000000000
#define MASK_SINGLE_NODE    0xC000000000000000


using namespace std;

namespace Workflow { namespace Cluster{

string MCL::get_description() {
	return "Markov clustering according to doi:10.1137/040608635";
}

void MCL::get_exp(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r){
#ifdef MCL_TIMINGS
	chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
#endif
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
#ifdef MCL_TIMINGS
	sparse_exp_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
#endif
}

void MCL::get_exp(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r){
	// TODO: at some r it may be more beneficial to diagnoalize in and only take the exponents of the eigenvalues
#ifdef MCL_TIMINGS
	chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
#endif
	if(r - (int) r == 0){
		*out = *in;
		for(uint32_t i=1; i<r; i++){ 
			(*out) *= (*in);
		}
#ifdef MCL_TIMINGS
		dense_int_exp_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
#endif
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
		if( thr > numeric_limits<float>::epsilon() ){
			out->noalias() = (V * d.real() * V.inverse()).real();
		}
#ifdef MCL_TIMINGS
		dense_gen_exp_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
#endif
	}
}

void MCL::get_gamma(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r){
#ifdef MCL_TIMINGS
	chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
#endif
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
	*out = out->pruned(1.0, numeric_limits<float>::epsilon());
#ifdef MCL_TIMINGS
	sparse_gamma_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
#endif
}

void MCL::get_gamma(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r){
#ifdef MCL_TIMINGS
	chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
#endif
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
#ifdef MCL_TIMINGS
	dense_gamma_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
#endif
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

void MCL::markov_process(Eigen::SparseMatrix<float>* m, float inflation, float expansion){
	for (uint32_t idiag=0; idiag<m->cols(); ++idiag){
		assert(abs(m->coeffRef(idiag, idiag)) > numeric_limits<float>::epsilon());
	}
	uint32_t iteration = 0;
	float diff_norm = numeric_limits<float>::max();
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
}

void MCL::markov_process(Eigen::MatrixXf* m, float inflation, float expansion){
	for (uint32_t idiag=0; idiag<m->cols(); ++idiag){
		assert(abs(m->coeffRef(idiag, idiag)) > numeric_limits<float>::epsilon());
	}
	uint32_t iteration = 0;
	float diff_norm = numeric_limits<float>::max();
	Eigen::MatrixXf msquared(m->rows(), m->cols());
	Eigen::MatrixXf m_update(m->rows(), m->cols());
	get_gamma(m, m, 1); // This is to get a matrix of random walks on the graph -> TODO: find out if something else is more suitable
	while( iteration < 100 && diff_norm > numeric_limits<float>::epsilon()*m->rows() ){
		get_exp(m, &msquared, expansion);
		get_gamma(&msquared, &m_update, inflation);
		*m -= m_update;
		diff_norm = m->norm();
		*m = m_update;
		iteration++;
	}
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
		format = "qcovhsp*scovhsp*pident";
	}
	config.output_format = {"clus", format};

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
	vector<uint32_t> sort_order(components.size());
	iota(sort_order.begin(), sort_order.end(), 0);
	sort(sort_order.begin(), sort_order.end(), [&](uint32_t i, uint32_t j){return indices[i].size() > indices[j].size();});

	uint32_t nComponents = count_if(indices.begin(), indices.end(), [](vector<uint32_t> v){ return v.size() > 0;});
	uint32_t nComponentsLt1 = count_if(indices.begin(), indices.end(), [](vector<uint32_t> v){ return v.size() > 1;});
	timer.finish();
	message_stream << "Found " << nComponentsLt1 <<" ("<<nComponents<< " incl. singletons) disconnected components" << endl;

	timer.go("Clustering components");
	bool full = false;
	// Note, we will access the clustering_result from several threads below and a vector does not guarantee thread-safety in these situations.
	// Note also, that the use of the disjoint_set structure guarantees that each thread will access a different part of the clustering_result
	uint64_t* clustering_result = new uint64_t[db->ref_header.sequences];
	fill(clustering_result, clustering_result+db->ref_header.sequences, 0UL);

	// note that the disconnected components are sorted by size
	const uint32_t chunk_size = 100;
	const uint32_t max_counter = nComponents;
	uint32_t nThreads = min(config.threads_, nComponents / chunk_size);
	const float inflation = (float) config.cluster_mcl_inflation;
	const float expansion = (float) config.cluster_mcl_expansion;

	// Collect some stats on the way
	float* max_sparsities = new float[nThreads];
	float* min_sparsities = new float[nThreads];
	atomic_uint n_clusters_found(0);
	atomic_uint component_counter(nThreads*chunk_size);
	atomic_uint n_dense_calculations(0);
	atomic_uint n_sparse_calculations(0);
	atomic_uint nClustersEq1(0);

	auto mcl_clustering = [&](const uint32_t iThr){
		uint32_t n_dense = 0;
		uint32_t n_sparse = 0;
		uint32_t n_singletons = 0;
		uint64_t cluster_id = iThr;
		uint32_t my_counter = iThr*chunk_size;
		float max_sparsity = 0.0f;
		float min_sparsity = 1.0f;
		while(my_counter < max_counter){
			for(uint32_t chunk_counter = my_counter; chunk_counter<min(my_counter+chunk_size, max_counter); chunk_counter++){
				uint32_t iComponent = sort_order[chunk_counter];
				vector<uint32_t> order = indices[iComponent];
				if(order.size() > 1){
					vector<Eigen::Triplet<float>> m = components[iComponent];
					assert(m.size() <= order.size() * order.size());
					float sparsity = 1.0-(1.0 * m.size()) / (order.size()*order.size());
					max_sparsity = max(max_sparsity, sparsity);
					min_sparsity = min(min_sparsity, sparsity);
					// map back to original ids
					vector<unordered_set<uint32_t>> list_of_sets;
					unordered_set<uint32_t> attractors;

					//TODO: a size limit for the dense matrix should control this as well
#ifdef MCL_TIMINGS
					chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
#endif
					if(sparsity >= config.cluster_mcl_sparsity_switch && expansion - (int) expansion == 0){ 
						n_sparse++;
						Eigen::SparseMatrix<float> m_sparse(order.size(), order.size());
						m_sparse.setFromTriplets(m.begin(), m.end());
#ifdef MCL_TIMINGS
						sparse_create_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
#endif
						markov_process(&m_sparse, inflation, expansion);
#ifdef MCL_TIMINGS
						chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
#endif
						LazyDisjointIntegralSet<uint32_t> disjointSet(m_sparse.cols());
						for (uint32_t k=0; k<m_sparse.outerSize(); ++k){
							for (Eigen::SparseMatrix<float>::InnerIterator it(m_sparse, k); it; ++it){
								assert(abs(it.value()) > numeric_limits<float>::epsilon());
								disjointSet.merge(it.row(), it.col());
								if(it.row() == it.col()){
									attractors.emplace(it.row());
								}
							}
						}
						list_of_sets = disjointSet.getListOfSets();
#ifdef MCL_TIMINGS
						sparse_list_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
#endif
					}
					else{
						n_dense++;
						Eigen::MatrixXf m_dense = Eigen::MatrixXf::Zero(order.size(), order.size());
						for(Eigen::Triplet<float> const & t : m){
							m_dense(t.row(), t.col()) = t.value();
						}
#ifdef MCL_TIMINGS
						dense_create_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
#endif
						markov_process(&m_dense, inflation, expansion);
#ifdef MCL_TIMINGS
						chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
#endif
						LazyDisjointIntegralSet<uint32_t> disjointSet(m_dense.cols());
						for (uint32_t icol=0; icol<m_dense.cols(); ++icol){
							for (uint32_t irow=0; irow<m_dense.rows(); ++irow){
								if( abs(m_dense(irow, icol)) > numeric_limits<float>::epsilon()){
									disjointSet.merge(irow, icol);
									if(irow == icol){
										attractors.emplace(irow);
									}
								}
							}
						}
						list_of_sets = disjointSet.getListOfSets();
#ifdef MCL_TIMINGS
						dense_list_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
#endif
					}
					for(unordered_set<uint32_t> subset : list_of_sets){
						assert(cluster_id < 0x3fffffffffffffff);
						for(uint32_t el : subset){
							const uint64_t mask = attractors.find(el) == attractors.end() ? MASK_NORMAL_NODE : MASK_ATTRACTOR_NODE;
							clustering_result[order[el]] = mask | cluster_id;
						}
						if(subset.size() == 1) n_singletons++;
						cluster_id += nThreads;
					}
				}
				else if (order.size() == 1){
					assert(cluster_id < 0x3fffffffffffffff);
					clustering_result[order[0]] = MASK_SINGLE_NODE | cluster_id;
					cluster_id += nThreads;
					n_singletons++;
				}
			}
			my_counter = component_counter.fetch_add(chunk_size, memory_order_relaxed);
		}
		n_clusters_found.fetch_add((cluster_id-iThr)/nThreads, memory_order_relaxed);
		n_dense_calculations.fetch_add(n_dense, memory_order_relaxed);
		n_sparse_calculations.fetch_add(n_sparse, memory_order_relaxed);
		nClustersEq1.fetch_add(n_singletons, memory_order_relaxed);
		max_sparsities[iThr] = max_sparsity;
		min_sparsities[iThr] = min_sparsity;
	};

	vector<thread> threads;
	for(uint32_t iThread = 0; iThread < nThreads ; iThread++ ) {
		min_sparsities[iThread] = 1.0f;
		max_sparsities[iThread] = 0.0f;
		threads.emplace_back(mcl_clustering, iThread);
	}

	float min_sparsity = 1.0f;
	float max_sparsity = 0.0f;
	for(uint32_t iThread = 0; iThread < nThreads ; iThread++ ) {
		threads[iThread].join();
		min_sparsity = min(min_sparsity, min_sparsities[iThread]);
		max_sparsity = max(max_sparsity, max_sparsities[iThread]);
	}
	delete[] min_sparsities;
	delete[] max_sparsities;

	timer.finish();
	message_stream << "Clusters found " << n_clusters_found - nClustersEq1 << " ("<< n_clusters_found <<" incl. singletons)" << endl; 
	message_stream << "\t max size " << indices[sort_order[0]].size() << endl;
	message_stream << "\t min sparsity " << min_sparsity << endl;
	message_stream << "\t max sparsity " << max_sparsity << endl;
	message_stream << "\t number of dense calculations " << n_dense_calculations << endl;
	message_stream << "\t number of sparse calculations " << n_sparse_calculations << endl;
#ifdef MCL_TIMINGS
	message_stream << "Time used for matrix creation: " 
		<< (sparse_create_time.load() + dense_create_time.load())/1000.0 
		<<" (sparse: " 
		<< sparse_create_time.load()/1000.0 
		<< ", dense: " 
		<< dense_create_time.load()/1000.0 <<")" << endl;
	message_stream << "Time used for exp: " 
		<< (sparse_exp_time.load() + dense_int_exp_time.load() + dense_gen_exp_time.load())/1000.0 
		<< " (sparse: " 
		<< sparse_exp_time.load()/1000.0 
		<< ", dense int: " 
		<< dense_int_exp_time.load()/1000.0 
		<< ", dense gen: " 
		<< dense_gen_exp_time.load()/1000.0 <<")" << endl; 
	message_stream << "Time used for gamma: " 
		<< (sparse_gamma_time.load() + dense_gamma_time.load())/1000.0 
		<< " (sparse: " 
		<< sparse_gamma_time.load()/1000.0 
		<< ", dense: " 
		<< dense_gamma_time.load()/1000.0 <<")" << endl;
	message_stream << "Time used for listing: " 
		<< (sparse_list_time.load() + dense_list_time.load())/1000.0 
		<< " (sparse: " 
		<< sparse_list_time.load()/1000.0 
		<< ", dense: " 
		<< dense_list_time.load()/1000.0 <<")" << endl;
#endif

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
		const uint64_t res = clustering_result[i];
		(*out) << blast_id(id) << '\t' << (res & (~MASK_INVERSE)) << '\t';
		switch(res & MASK_INVERSE){
			case MASK_SINGLE_NODE:
				(*out) << 's';
				break;
			case MASK_ATTRACTOR_NODE:
				(*out) << 'a';
				break;
			case MASK_NORMAL_NODE:
				(*out) << 'n';
				break;
			default:
				// Note that this is a sanity check that we have touched all elements in clustering_result
				(*out) << 'u' << "\tERROR: This should not occur";
		}
		(*out) << endl;
		id.clear();
		seq.clear();
	}
	if (out != &cout) delete out;
	db->close();
	delete[] clustering_result;
	timer.finish();
}
}}
