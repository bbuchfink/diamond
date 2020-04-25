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

#pragma once
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
#include "../run/workflow.h"
#include "disjoint_set.h"
#include "../util/io/consumer.h"
#include "../util/algo/algo.h"
#include "../basic/statistics.h"
#include "../util/log_stream.h"
#include "../dp/dp.h"
#include "cluster.h"

using namespace std;

namespace Workflow { namespace Cluster{
class MCL: public ClusteringAlgorithm {
private: 
	Eigen::SparseMatrix<float> get_gamma(Eigen::SparseMatrix<float>* m, float r);
	vector<unordered_set<uint32_t>> get_list(Eigen::SparseMatrix<float>* m);
	vector<unordered_set<uint32_t>> markov_process(Eigen::SparseMatrix<float>* m);
	void run_clustering();
public:
	void run(){
		run_clustering();
	}
	static string get_key(){
		return "mcl";
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
		return make_pair(move(indices), move(components));
	}
	Eigen::SparseMatrix<T> getMatrix(){
		Eigen::SparseMatrix<T> m(n,n);
		m.setFromTriplets(data.begin(), data.end());
		return m;
	}
};
}}
