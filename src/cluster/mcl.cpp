#include "mcl.h"

using namespace std;

namespace Workflow { namespace Cluster{

Eigen::SparseMatrix<float> MCL::get_gamma(Eigen::SparseMatrix<float>* m, float r){
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

vector<unordered_set<uint32_t>> MCL::get_list(Eigen::SparseMatrix<float>* m){
	LazyDisjointIntegralSet<uint32_t> disjointSet(m->cols());
	for (uint32_t k=0; k<m->outerSize(); ++k){
		for (Eigen::SparseMatrix<float>::InnerIterator it(*m, k); it; ++it){
			disjointSet.merge(it.row(), it.col());
		}
	}
	return disjointSet.getListOfSets();
}

vector<unordered_set<uint32_t>> MCL::markov_process(Eigen::SparseMatrix<float>* m){
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

void MCL::run_clustering() {
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
	config.no_self_hits = true;
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
	// 
	db->close();
}
}}
