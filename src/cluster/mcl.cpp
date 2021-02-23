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
#include <iomanip>
#include "mcl.h"
#include "sparse_matrix_stream.h"

#define MASK_INVERSE        0xC000000000000000
#define MASK_NORMAL_NODE    0x4000000000000000
#define MASK_ATTRACTOR_NODE 0x8000000000000000
#define MASK_SINGLE_NODE    0xC000000000000000


using namespace std;

namespace Workflow { namespace Cluster{

string MCL::get_description() {
	return "Markov clustering according to doi:10.1137/040608635";
}

inline uint64_t get_idx(uint32_t i_idx, uint32_t n_cols, uint32_t j_idx){
	return i_idx * n_cols + j_idx;
}

inline void get_indices (uint64_t combined, uint32_t* i, uint32_t n_cols, uint32_t* j){
	*j = combined % n_cols;
	*i = ( combined - *j) / n_cols;
}

vector<Eigen::Triplet<float>> multiply(Eigen::SparseMatrix<float>* a, Eigen::SparseMatrix<float>* b , uint32_t iThr, uint32_t nThr){
	uint32_t n_cols = b->cols();
	uint32_t n_rows = a->rows();
	std::vector<float> result_col(n_rows, 0.0);
	vector<Eigen::Triplet<float>> data;
	for (uint32_t j=0; j<n_cols; ++j) {
		if(j%nThr == iThr){
			fill(result_col.begin(), result_col.end(), 0.0);
			for (Eigen::SparseMatrix<float>::InnerIterator  rhsIt(*b, j); rhsIt; ++rhsIt) {
				const float y = rhsIt.value();
				const uint32_t k = rhsIt.row();
				for (Eigen::SparseMatrix<float>::InnerIterator lhsIt(*a, k); lhsIt; ++lhsIt) {
					const uint32_t i = lhsIt.row();
					const float x = lhsIt.value();
					result_col[i] += x*y;
				}
			}
			for (uint32_t i=0; i<n_rows; ++i) {
				if(result_col[i] > numeric_limits<float>::epsilon()) {
					data.emplace_back(i, j, result_col[i]);
				}
			}
		}
	}
	return data;
}

void MCL::get_exp(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r, uint32_t nThr){
	chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
	if( r - (int) r == 0 ){
		// TODO: at some r it may be more beneficial to diagnoalize in and only take the exponents of the eigenvalues
		*out = *in;
		for(uint32_t i=1; i<r; i++){ 
			// Note: this is a workaround for a crash in Eigen
			vector<Eigen::Triplet<float>> data;
			std::mutex m;
			auto mult = [&](const uint32_t iThr){
				vector<Eigen::Triplet<float>> t_data = multiply(in, out, iThr, nThr);
				m.lock();
				data.insert(data.end(), t_data.begin(), t_data.end());
				m.unlock();
			};

			vector<thread> threads;
			for(uint32_t iThread = 0; iThread < nThr ; iThread++) {
				threads.emplace_back(mult, iThread);
			}

			for(uint32_t iThread = 0; iThread < nThr ; iThread++) {
				threads[iThread].join();
			}
			out->setZero();
			out->setFromTriplets(data.begin(), data.end(), [] (const float&, const float &b) { return b; });
			*out = out->pruned(1.0, numeric_limits<float>::epsilon());
			//(*out) = (*out)*(*in);
		}
	}
	else{ 
		throw runtime_error(" Eigen does not provide an eigenvalue solver for sparse matrices");
	}
	*out = out->pruned();
	sparse_exp_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
}

void MCL::get_exp(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r){
	// TODO: at some r it may be more beneficial to diagnoalize in and only take the exponents of the eigenvalues
	chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
	if(r - (int) r == 0){
		*out = *in;
		for(uint32_t i=1; i<r; i++){ 
			(*out) *= (*in);
		}
		//*out = in->pow(r)
		dense_int_exp_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
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
		dense_gen_exp_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
	}
}


vector<Eigen::Triplet<float>> gammaIze(Eigen::SparseMatrix<float>* in, float r, uint32_t iThr, uint32_t nThr){
	vector<Eigen::Triplet<float>> data;
	for (uint32_t k=0; k<in->outerSize(); ++k){
		if( k%nThr == iThr){
			float colSum = 0.0f;
			for (Eigen::SparseMatrix<float>::InnerIterator it(*in, k); it; ++it){
				colSum += pow(it.value(),r);
			}
			for (Eigen::SparseMatrix<float>::InnerIterator it(*in, k); it; ++it){
				float val = pow(it.value(), r) / colSum;
				if(val > numeric_limits<float>::epsilon()) {
					data.emplace_back(it.row(), it.col(), val);
				}
			}
		}
	}
	return data;
}

void MCL::get_gamma(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r, uint32_t nThr){
	chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
	vector<Eigen::Triplet<float>> data;
	std::mutex m;
	auto mult = [&](const uint32_t iThr){
		vector<Eigen::Triplet<float>> t_data = gammaIze(in, r, iThr, nThr);
		m.lock();
		data.insert(data.end(), t_data.begin(), t_data.end());
		m.unlock();
	};

	vector<thread> threads;
	for(uint32_t iThread = 0; iThread < nThr ; iThread++) {
		threads.emplace_back(mult, iThread);
	}

	for(uint32_t iThread = 0; iThread < nThr ; iThread++) {
		threads[iThread].join();
	}
	out->setZero();
	out->setFromTriplets(data.begin(), data.end(), [] (const float&, const float &b) { return b; });
	*out = out->pruned(1.0, numeric_limits<float>::epsilon());
	sparse_gamma_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
}

float get_norm(Eigen::SparseMatrix<float>* in, uint32_t nThr){
	vector<float> data(nThr);
	fill(data.begin(), data.end(), 0.0);
	auto norm = [&](const uint32_t iThr){
		for (uint32_t k=0; k<in->outerSize(); ++k){
			if( k%nThr == iThr){
				for (Eigen::SparseMatrix<float>::InnerIterator it(*in, k); it; ++it){
					data[iThr] += pow(it.value(),2.0);
				}
			}
		}
	};

	vector<thread> threads;
	for(uint32_t iThread = 0; iThread < nThr ; iThread++) {
		threads.emplace_back(norm, iThread);
	}

	for(uint32_t iThread = 0; iThread < nThr ; iThread++) {
		threads[iThread].join();
	}
	return pow(accumulate(data.begin(), data.end(), 0.0f), 0.5f);
}

void MCL::get_gamma(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r){
	chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
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
	dense_gamma_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
}

void MCL::markov_process(Eigen::SparseMatrix<float>* m, float inflation, float expansion, uint32_t max_iter, function<uint32_t()> getThreads){
	for (uint32_t idiag=0; idiag<m->cols(); ++idiag){
		assert(abs(m->coeffRef(idiag, idiag)) > numeric_limits<float>::epsilon());
	}
	uint32_t iteration = 0;
	float diff_norm = numeric_limits<float>::max();
	Eigen::SparseMatrix<float> msquared(m->rows(), m->cols());
	Eigen::SparseMatrix<float> m_update(m->rows(), m->cols());
	get_gamma(m, m, 1, getThreads()); // This is to get a matrix of random walks on the graph -> TODO: find out if something else is more suitable
	while( iteration < max_iter && diff_norm > 1e-6*m->rows() ){
		get_exp(m, &msquared, expansion, getThreads());
		get_gamma(&msquared, &m_update, inflation, getThreads());
		*m -= m_update;
		diff_norm = get_norm(m, getThreads());
		*m = m_update;
		iteration++;
	}
	if( iteration == max_iter ){
		failed_to_converge++;
	}
}

void MCL::markov_process(Eigen::MatrixXf* m, float inflation, float expansion, uint32_t max_iter){
	for (uint32_t idiag=0; idiag<m->cols(); ++idiag){
		assert(abs(m->coeffRef(idiag, idiag)) > numeric_limits<float>::epsilon());
	}
	uint32_t iteration = 0;
	float diff_norm = numeric_limits<float>::max();
	Eigen::MatrixXf msquared(m->rows(), m->cols());
	Eigen::MatrixXf m_update(m->rows(), m->cols());
	get_gamma(m, m, 1); // This is to get a matrix of random walks on the graph -> TODO: find out if something else is more suitable
	while( iteration < max_iter && diff_norm > numeric_limits<float>::epsilon()*m->rows() ){
		get_exp(m, &msquared, expansion);
		get_gamma(&msquared, &m_update, inflation);
		*m -= m_update;
		diff_norm = m->norm();
		*m = m_update;
		iteration++;
	}
	if( iteration == max_iter ){
		failed_to_converge++;
	}
}

void MCL::print_stats(uint64_t nElements, uint32_t nComponents, uint32_t nComponentsLt1, vector<uint32_t>& sort_order, vector<vector<uint32_t>>& indices, vector<vector<Eigen::Triplet<float>>>& components){
	const uint32_t chunk_size = config.cluster_mcl_chunk_size;
	const uint32_t nThreads = min(config.threads_, nComponents / chunk_size);
	const float inflation = (float) config.cluster_mcl_inflation;
	const float expansion = (float) config.cluster_mcl_expansion;
	vector<float> sparsities(nComponentsLt1);
	uint64_t icomp = 0;
	generate(sparsities.begin(), sparsities.end(), [&icomp, &sort_order, &indices, &components](){ 
				const uint32_t iComponent = sort_order[icomp++];
				const uint32_t ssize = components[iComponent].size();
				const uint64_t dsize = indices[iComponent].size() * indices[iComponent].size();
				return 1.0-(1.0 * ssize) / dsize;
			});
	vector<float> neighbors(nComponentsLt1);
	icomp = 0;
	generate(neighbors.begin(), neighbors.end(), [&icomp, &sort_order, &indices, &components](){ 
				const uint32_t iComponent = sort_order[icomp++];
				uint32_t neighbors = 0;
				for(auto t: components[iComponent]){
					if(t.row() != t.col()) neighbors++;
				}
				return neighbors / ((float) indices.size());
			});

	// rough memory calc before sorting sparsities
	vector<float> memories(nComponentsLt1);
	icomp = 0;
	generate(memories.begin(), memories.end(), [&](){ 
				const uint32_t iComponent = sort_order[icomp++];
				if(sparsities[icomp - 1] >= config.cluster_mcl_sparsity_switch){ 
					return indices[iComponent].size() * (1+pow(neighbors[icomp], expansion)) * (2 * sizeof(uint32_t) + sizeof(float));
				}
				else{
					return (float) sizeof(float) * indices[iComponent].size() * indices[iComponent].size();
				}
			});
	sort(memories.rbegin(), memories.rend());
	float mem_req = nElements * (2 * sizeof(uint32_t) + sizeof(float));
	for(uint32_t iThr = 0; iThr < nThreads; iThr++){
		mem_req += 3*memories[iThr]; // we need three matrices of this size
	}
	
	float neighbors_of_max = neighbors[0];
	float neighbors_of_med = nComponentsLt1 > 0 && nComponentsLt1 % 2 == 0 ? 
		(neighbors[nComponentsLt1/2] + neighbors[nComponentsLt1/2+1]) / 2.0 : 
		neighbors[nComponentsLt1/2+1] ;
	float neighbors_of_min = neighbors[nComponentsLt1-1];

	float sparsity_of_max = sparsities[0];
	float sparsity_of_med = nComponentsLt1 > 0 && nComponentsLt1 % 2 == 0 ? 
		(sparsities[nComponentsLt1/2] + sparsities[nComponentsLt1/2+1]) / 2.0 : 
		sparsities[nComponentsLt1/2+1] ;
	float sparsity_of_min = sparsities[nComponentsLt1-1];

	sort(sparsities.rbegin(), sparsities.rend());
	sort(neighbors.rbegin(), neighbors.rend());

	float med_sparsity = nComponentsLt1 > 0 && nComponentsLt1 % 2 == 0 ? 
		(sparsities[nComponentsLt1/2] + sparsities[nComponentsLt1/2+1]) / 2.0 : 
		sparsities[nComponentsLt1/2+1] ;
	float med_neighbors = nComponentsLt1 > 0 && nComponentsLt1 % 2 == 0 ? 
		(neighbors[nComponentsLt1/2] + neighbors[nComponentsLt1/2+1]) / 2.0 : 
		neighbors[nComponentsLt1/2+1] ;
	float median1 = nComponents > 0 && nComponents % 2 == 0 ? 
		(indices[sort_order[nComponents/2]].size() + indices[sort_order[nComponents/2+1]].size()) / 2.0 : 
		indices[sort_order[nComponents/2+1]].size();
	float median2 = nComponentsLt1 > 0 && nComponentsLt1 % 2 == 0 ? 
		(indices[sort_order[nComponentsLt1/2]].size() + indices[sort_order[nComponentsLt1/2+1]].size()) / 2.0 : 
		indices[sort_order[nComponentsLt1/2+1]].size();

	message_stream << "Number of DIAMOND hits:          " << nElements << endl;
	message_stream << "Number of independet components: " << nComponentsLt1 <<" ("<<nComponents<< " incl. singletons)" << endl;
	message_stream << "Component size information: "<< endl;
	message_stream << "\tmax. : " << setw(12) << indices[sort_order[0]].size() << endl;
	message_stream << "\tmed. : " << setw(12) << median2 << " ("<<median1<< " incl. singletons)"<<endl;
	message_stream << "\tmin. : " << setw(12) << indices[sort_order[nComponentsLt1-1]].size() << " ("<<indices[sort_order[nComponents-1]].size()<< " incl. singletons)"<<endl;
	message_stream << "Sparsity of components (excluding singletons): "<< endl;
	message_stream << "\tmax. : " << setw(12) << sparsities[0]                << " - at max. size: " << setw(12) << sparsity_of_max << endl;
	message_stream << "\tmed. : " << setw(12) << med_sparsity                 << " - at med. size: " << setw(12) << sparsity_of_med << endl;
	message_stream << "\tmin. : " << setw(12) << sparsities[nComponentsLt1-1] << " - at min. size: " << setw(12) << sparsity_of_min << endl;
	message_stream << "Average number of neighbors in components (excluding singletons): "<< endl;
	message_stream << "\tmax. : " << setw(12) << neighbors[0]                << " - at max. size: " << setw(12) << neighbors_of_max << endl;
	message_stream << "\tmed. : " << setw(12) << med_neighbors               << " - at med. size: " << setw(12) << neighbors_of_med << endl;
	message_stream << "\tmin. : " << setw(12) << neighbors[nComponentsLt1-1] << " - at min. size: " << setw(12) << neighbors_of_min << endl;
	message_stream << "Rough memory requirements: ";
	string uints[] = {"", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};
	uint32_t idx = ((uint32_t)log(mem_req)/log(1024));
	message_stream << setprecision(2) << mem_req/pow(1024,idx) << uints[idx] << endl;
}

SparseMatrixStream<float>* get_graph_handle(SequenceFile& db){
	if(config.cluster_mcl_restart){
		task_timer timer;
		timer.go("Reading cluster checkpoint file");
		SparseMatrixStream<float>* ms = SparseMatrixStream<float>::fromFile(config.cluster_mcl_graph_file, config.chunk_size);
		timer.finish();
		ms->done();
		return ms;
	}
	config.command = Config::blastp;
	config.no_self_hits = false;
	string format = config.cluster_similarity;
	if(format.empty()){
		format = "qcovhsp*scovhsp*pident";
	}
	config.output_format = {"clus", format};
	Workflow::Search::Options opt;
	opt.db = &db;
	opt.self = true;
	SparseMatrixStream<float>* ms = new SparseMatrixStream<float>(db.sequence_count(), config.cluster_mcl_graph_file);
	if(config.chunk_size > 0){
		ms->set_max_mem(config.chunk_size);
	}
	opt.consumer = ms;
	Workflow::Search::run(opt);
	ms->done();
	return ms;
}

void MCL::run(){
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
	// unordered_set<uint32_t> s ( {183, 239, 23834, 3399, 484, 2} );
	// LazyDisjointTypeSet<uint32_t> ts(&s);

	// ts.merge(183, 3399);
	// ts.merge(23834, 2);
	// ts.merge(183, 484);
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

	if (config.database == "") throw runtime_error("Missing parameter: database file (--db/-d)");
	unique_ptr<SequenceFile> db(SequenceFile::auto_create());
	statistics.reset();
	SparseMatrixStream<float>* ms = get_graph_handle(*db);
	task_timer timer;
	timer.go("Computing independent components");
	vector<vector<uint32_t>> indices = ms->get_indices();
	ms->clear_disjoint_set();
	uint64_t nElements = ms->getNumberOfElements();
	uint32_t nComponents = count_if(indices.begin(), indices.end(), [](vector<uint32_t> v){ return v.size() > 0;});
	uint32_t nComponentsLt1 = count_if(indices.begin(), indices.end(), [](vector<uint32_t> v){ return v.size() > 1;});

	// Sort to get the largest chunks first
	vector<uint32_t> sort_order(indices.size());
	iota(sort_order.begin(), sort_order.end(), 0);
	sort(sort_order.begin(), sort_order.end(), [&](uint32_t i, uint32_t j){return indices[i].size() > indices[j].size();});
	timer.finish();
	//print_stats(nElements, nComponents, nComponentsLt1, sort_order, indices, components);

	timer.go("Clustering components");
	// Note, we will access the clustering_result from several threads below and a vector does not guarantee thread-safety in these situations.
	// Note also, that the use of the disjoint_set structure guarantees that each thread will access a different part of the clustering_result
	uint64_t* clustering_result = new uint64_t[db->sequence_count()];
	fill(clustering_result, clustering_result+db->sequence_count(), 0UL);

	// note that the disconnected components are sorted by size
	const uint32_t chunk_size = config.cluster_mcl_chunk_size;
	const uint32_t nThreads = min(config.threads_, nComponents / chunk_size);
	ms->allocate_read_buffer(nThreads);
	const float inflation = (float) config.cluster_mcl_inflation;
	const float expansion = (float) config.cluster_mcl_expansion;
	const uint32_t max_counter = nComponents;
	const uint32_t max_iter = config.cluster_mcl_max_iter;
	uint32_t max_job_size = 0;
	for(uint32_t i=0; i<chunk_size; i++) max_job_size+=indices[sort_order[i]].size();
	
	// Collect some stats on the way
	uint32_t* jobs_per_thread = new uint32_t[nThreads];
	float* time_per_thread = new float[nThreads];
	atomic_uint n_clusters_found(0);
	atomic_uint component_counter(nThreads*chunk_size);
	atomic_uint n_dense_calculations(0);
	atomic_uint n_sparse_calculations(0);
	atomic_uint nClustersEq1(0);
	atomic_uint threads_done(0);

	auto mcl_clustering = [&](const uint32_t iThr){
		chrono::high_resolution_clock::time_point thread_start = chrono::high_resolution_clock::now();
		uint32_t n_dense = 0;
		uint32_t n_sparse = 0;
		uint32_t n_singletons = 0;
		uint32_t n_jobs_done = 0;
		uint64_t cluster_id = iThr;
		uint32_t my_counter = iThr*chunk_size;
		uint32_t my_chunk_size = chunk_size;
		while(my_counter < max_counter){
			unordered_set<uint32_t> all;
			uint32_t upper_limit =  min(my_counter+my_chunk_size, max_counter);
			vector<vector<uint32_t>*> loc_i;
			my_chunk_size = 0;
			for(uint32_t chunk_counter = my_counter; chunk_counter<upper_limit; chunk_counter++){
				loc_i.push_back(&indices[sort_order[chunk_counter]]);
				my_chunk_size += indices[sort_order[chunk_counter]].size();
			}
			vector<vector<Eigen::Triplet<float>>> loc_c = ms->collect_components(&loc_i, iThr);
			uint32_t ichunk = 0;
			for(uint32_t chunk_counter = my_counter; chunk_counter<upper_limit; chunk_counter++){
				n_jobs_done++;
				vector<uint32_t>* order = loc_i[ichunk];
				if(order->size() > 1){
					vector<Eigen::Triplet<float>>* m = &loc_c[ichunk];
					assert(m->size() <= order->size() * order->size());
					// map back to original ids
					vector<unordered_set<uint32_t>> list_of_sets;
					unordered_set<uint32_t> attractors;

					const float sparsity = 1.0-(1.0 * m->size()) / (order->size()*order->size());
					chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
					//TODO: a size limit for the dense matrix should control this as well
					if(sparsity >= config.cluster_mcl_sparsity_switch && expansion - (int) expansion == 0){ 
						n_sparse++;
						Eigen::SparseMatrix<float> m_sparse(order->size(), order->size());
						m_sparse.setFromTriplets(m->begin(), m->end());
						sparse_create_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
						auto getThreads = [&threads_done, &nThreads, &iThr](){
							uint32_t td=threads_done.load();
							uint32_t rt=nThreads-td;
							return 1 + td / rt + (iThr < (td % rt) ? 1 : 0);
						};
						markov_process(&m_sparse, inflation, expansion, max_iter, getThreads);
						chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
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
						sparse_list_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
					}
					else{
						n_dense++;
						Eigen::MatrixXf m_dense = Eigen::MatrixXf::Zero(order->size(), order->size());
						for(Eigen::Triplet<float> const & t : *m){
							m_dense(t.row(), t.col()) = t.value();
						}
						dense_create_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
						markov_process(&m_dense, inflation, expansion, max_iter);
						chrono::high_resolution_clock::time_point t = chrono::high_resolution_clock::now();
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
						dense_list_time += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
					}
					for(unordered_set<uint32_t> subset : list_of_sets){
						assert(cluster_id < 0x3fffffffffffffff);
						for(uint32_t el : subset){
							const uint64_t mask = attractors.find(el) == attractors.end() ? MASK_NORMAL_NODE : MASK_ATTRACTOR_NODE;
							clustering_result[(*order)[el]] = mask | cluster_id;
						}
						if(subset.size() == 1) n_singletons++;
						cluster_id += nThreads;
					}
				}
				else if (order->size() == 1){
					assert(cluster_id < 0x3fffffffffffffff);
					clustering_result[(*order)[0]] = MASK_SINGLE_NODE | cluster_id;
					cluster_id += nThreads;
					n_singletons++;
				}
				ichunk++;
			}
			my_chunk_size = 0;
			for(auto const i : loc_i) my_chunk_size += i->size();
			my_chunk_size = loc_i.size() * (max_job_size / my_chunk_size);
			my_counter = component_counter.fetch_add(my_chunk_size, memory_order_relaxed);
		}
		n_clusters_found.fetch_add((cluster_id-iThr)/nThreads, memory_order_relaxed);
		n_dense_calculations.fetch_add(n_dense, memory_order_relaxed);
		n_sparse_calculations.fetch_add(n_sparse, memory_order_relaxed);
		nClustersEq1.fetch_add(n_singletons, memory_order_relaxed);
		jobs_per_thread[iThr] = n_jobs_done;
		threads_done.fetch_add(1, memory_order_relaxed);
		time_per_thread[iThr] = (chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - thread_start).count()) / 1000.0;
	};

	vector<thread> threads;
	for(uint32_t iThread = 0; iThread < nThreads ; iThread++) {
		threads.emplace_back(mcl_clustering, iThread);
	}

	for(uint32_t iThread = 0; iThread < nThreads ; iThread++) {
		threads[iThread].join();
	}
	ms->release_read_buffer();
	delete ms;
	timer.finish();

	// Output stats
	message_stream << "Jobs per thread: ";
	for(uint32_t iThread = 0; iThread < nThreads ; iThread++) {
		message_stream << " " << setw(8) << right << jobs_per_thread[iThread];
	}
	message_stream << endl; 
	message_stream << "Time per thread: ";
	for(uint32_t iThread = 0; iThread < nThreads ; iThread++) {
		message_stream << " " << setw(8) << right << time_per_thread[iThread];
	}
	message_stream << endl;
	delete[] time_per_thread;
	delete[] jobs_per_thread;

	message_stream << "Clusters found " << n_clusters_found - nClustersEq1 << " ("<< n_clusters_found <<" incl. singletons)" << endl; 
	message_stream << "\t number of failed calculations " << failed_to_converge.load() << endl;
	message_stream << "\t number of dense calculations " << n_dense_calculations << endl;
	message_stream << "\t number of sparse calculations " << n_sparse_calculations << endl;
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

	timer.go("Cluster output");
	ostream *out = config.output_file.empty() ? &cout : new ofstream(config.output_file.c_str());
	vector<Letter> seq;
	string id;
	db->init_seq_access();
	Hsp hsp;
	size_t n;
	out->precision(3);
	for (int i = 0; i < (int)db->sequence_count(); ++i) {
		db->read_seq(seq, id);
		const uint64_t cid = (~MASK_INVERSE) & clustering_result[i];
		const uint64_t lab = MASK_INVERSE & clustering_result[i];
		(*out) << blast_id(id) << '\t' ;
		switch(lab){
			case MASK_SINGLE_NODE:
				(*out) << cid << '\t'<< 's';
				break;
			case MASK_ATTRACTOR_NODE:
				(*out) << cid << '\t'<< 'a';
				break;
			case MASK_NORMAL_NODE:
				(*out) << cid << '\t'<< 'n';
				break;
			default:
				// Note that this is a sanity check that we have touched all elements in clustering_result
				(*out) << -1 << '\t' << 'u';
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
