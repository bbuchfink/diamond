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

#pragma once
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <stdio.h>
#include <memory>
#include <fstream>
#include <iostream>
#include <limits>
#include <atomic>
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
	void print_stats(uint64_t nElements, uint32_t nComponents, uint32_t nComponentsLt1, vector<uint32_t>& sort_order, vector<vector<uint32_t>>& indices, vector<vector<Eigen::Triplet<float>>>& components);
	void get_exp(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r);
	void get_exp(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r);
	void get_gamma(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r);
	void get_gamma(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r);
	void markov_process(Eigen::SparseMatrix<float>* m, float inflation, float expansion, uint32_t max_iter);
	void markov_process(Eigen::MatrixXf* m, float inflation, float expansion, uint32_t max_iter);
	atomic_ullong failed_to_converge = {0};
	atomic_ullong sparse_create_time = {0};
	atomic_ullong dense_create_time = {0};
	atomic_ullong sparse_exp_time = {0};
	atomic_ullong dense_int_exp_time = {0};
	atomic_ullong dense_gen_exp_time = {0};
	atomic_ullong sparse_gamma_time = {0};
	atomic_ullong dense_gamma_time = {0};
	atomic_ullong sparse_list_time = {0};
	atomic_ullong dense_list_time = {0};
public:
	~MCL(){};
	void run();
	string get_description();
	static string get_key(){
		return "mcl";
	}
};

template <typename T>
class SparseMatrixStream : public Consumer {
	struct CoordinateCmp {
		bool operator()(const Eigen::Triplet<T>& lhs, const Eigen::Triplet<T>& rhs) const { 
			return lhs.row() == rhs.row() ? lhs.col() < rhs.col() : lhs.row() < rhs.row() ; 
		}
	};
	size_t n;
	float max_size;
	set<Eigen::Triplet<T>, CoordinateCmp> data;
	LazyDisjointSet<uint32_t>* disjointSet;
	ofstream* os;
	inline void write_triplet(const uint32_t* const query, const uint32_t* const subject, const T* const value){
		if(os){
			os->write((char*) query, sizeof(uint32_t));
			os->write((char*) subject, sizeof(uint32_t));
			double val = *value;
			os->write((char*) &val, sizeof(double));
		}
	}
	void dump(){
		if(os){
			auto c = getComponents();
			vector<vector<uint32_t>> indices = get<0>(c);
			vector<vector<Eigen::Triplet<T>>> components = get<1>(c);
			for(uint32_t iComponent = 0; iComponent<indices.size(); iComponent++){
				if( components[iComponent].size() > 0 ){
					const uint32_t first_index = indices[iComponent][0];
					const uint32_t size = components[iComponent].size();
					os->write((char*) &first_index, sizeof(uint32_t));
					os->write((char*) &size, sizeof(uint32_t));
					for(auto& t : components[iComponent]){
						const uint32_t irow = indices[iComponent][t.row()];
						const uint32_t icol = indices[iComponent][t.col()];
						const T val = t.value();
						write_triplet(&irow, &icol, &val);
					}
				}
			}
		}
	}
	virtual void consume(const char *ptr, size_t n) override {
		const char *end = ptr + n;
		if (n >= sizeof(uint32_t)) {
			while (ptr  < end) {
				const uint32_t query = *(uint32_t*) ptr;
				ptr += sizeof(uint32_t);
				const uint32_t subject = *(uint32_t*) ptr;
				ptr += sizeof(uint32_t);
				const double value = *(double*) ptr;
				ptr += sizeof(double);
				Eigen::Triplet<T> t(query, subject, value);
				auto it = data.find(t);
				if(it == data.end()){
					data.emplace(move(t));
					disjointSet->merge(query, subject);
				}
				else if (t.value() > it->value()){
					data.emplace_hint(data.erase(it), move(t));
				}
				if(os && data.size()*(2*sizeof(uint32_t)+sizeof(T))/(1024ULL*1024*1024) >= max_size){
					dump();
					data.clear();
				}
			}
		}
	}

public:
	SparseMatrixStream(size_t n, string graph_file_name){
		this->n = n;
		disjointSet = new LazyDisjointIntegralSet<uint32_t>(n); 
		this->max_size = 5.0; // 5 GB default flush
		if(graph_file_name.empty()){
			os = nullptr;
		}
		else{
			os = new ofstream(graph_file_name, ios::out | ios::binary);
			os->write((char*) &n, sizeof(size_t));
			const uint32_t indexversion = 0;
			os->write((char*) &indexversion, sizeof(uint32_t));
		}
	}
	SparseMatrixStream(size_t n){
		this->n = n;
		disjointSet = new LazyDisjointIntegralSet<uint32_t>(n); 
		os = nullptr;
	}

	~SparseMatrixStream(){
		dump();
		data.clear();
		if(disjointSet) delete disjointSet;
		if(os){
			os->close();
			delete os;
		}
	}

	void set_max_mem(float max_size){
		this->max_size = max_size;
	}

	static SparseMatrixStream fromFile(string graph_file_name){
		ifstream in(graph_file_name, ios::in | ios::binary);
		if(!in){
			throw runtime_error("Cannot read the graph file");
		}
		size_t n;
		in.read((char*) &n, sizeof(size_t));
		uint32_t indexversion;
		in.read((char*) &indexversion, sizeof(uint32_t));
		if(indexversion != 0){
			throw runtime_error("This file cannot be read");
		}
		SparseMatrixStream sms(n);
		size_t unit_size = 2*sizeof(uint32_t)+sizeof(double);
		unsigned long long buffer_size = (2ULL*1024*1024*1024 / unit_size) * unit_size;
		char* buffer = new char[buffer_size];
		while(in.good()){
			uint32_t first_component;
			in.read((char*) &first_component, sizeof(uint32_t));
			uint32_t size;
			in.read((char*) &size, sizeof(uint32_t));
			if(size*unit_size > buffer_size){
				delete[] buffer;
				buffer_size = size*unit_size;
				buffer = new char[buffer_size];
			}
			in.read(buffer, size*unit_size);
			sms.consume(buffer, size*unit_size);
		}
		in.close();
		return sms;
	}

	static pair<vector<vector<uint32_t>>, vector<vector<Eigen::Triplet<T>>>> collect_components(unordered_set<uint32_t>& set, string graph_file_name){
		ifstream in(graph_file_name, ios::in | ios::binary);
		if(!in){
			throw runtime_error("Cannot read the graph file");
		}
		size_t n;
		in.read((char*) &n, sizeof(size_t));
		uint32_t indexversion;
		in.read((char*) &indexversion, sizeof(uint32_t));
		if(indexversion != 0){
			throw runtime_error("This file cannot be read");
		}
		SparseMatrixStream sms(n);
		size_t unit_size = 2*sizeof(uint32_t)+sizeof(double);
		unsigned long long buffer_size = (2ULL*1024*1024*1024 / unit_size) * unit_size;
		char* buffer = new char[buffer_size];
		while(in.good()){
			uint32_t first_component;
			in.read((char*) &first_component, sizeof(uint32_t));
			uint32_t size;
			in.read((char*) &size, sizeof(uint32_t));
			if(size*unit_size > buffer_size){
				delete[] buffer;
				buffer_size = size*unit_size;
				buffer = new char[buffer_size];
			}
			if(set.count(first_component) == 1){
				in.read(buffer, size*unit_size);
				sms.consume(buffer, size*unit_size);
			}
			else{
				in.ignore(size*unit_size);
			}
		}
		in.close();
		return sms.getComponents();
	}

	pair<vector<vector<uint32_t>>, vector<vector<Eigen::Triplet<T>>>> getComponents(){
		vector<unordered_set<uint32_t>> sets = disjointSet->getListOfSets();
		vector<uint32_t> indexToSetId(this->n, sets.size());
		for(uint32_t iset = 0; iset < sets.size(); iset++){
			for(uint32_t index : sets[iset]){
				indexToSetId[index] = iset;
			}
		}
		vector<vector<Eigen::Triplet<T>>> split(sets.size());
		for(Eigen::Triplet<T> t : data){
			uint32_t iset = indexToSetId[t.row()];
			assert( iset == indexToSetId[t.col()]);
			split[iset].emplace_back(t.row(), t.col(), t.value());
		}
		vector<vector<uint32_t>> indices;
		vector<vector<Eigen::Triplet<T>>> components;
		for(uint32_t iset = 0; iset < sets.size(); iset++){
			if(split[iset].size() > 0){
				vector<uint32_t> order(sets[iset].begin(), sets[iset].end());
				map<uint32_t, uint32_t> index_map;
				uint32_t iel = 0;
				for(uint32_t const & el: order){
					index_map.emplace(el, iel++);
				}
				vector<Eigen::Triplet<T>> remapped;
				for(Eigen::Triplet<T> const & t : split[iset]){
					remapped.emplace_back(index_map[t.row()], index_map[t.col()], t.value());
				}
				remapped.shrink_to_fit();
				components.emplace_back(move(remapped));
				indices.emplace_back(move(order));
			}
		}
		return make_pair(move(indices), move(components));
	}

	Eigen::SparseMatrix<T> getMatrix(){
		Eigen::SparseMatrix<T> m(n,n);
		m.setFromTriplets(data.begin(), data.end());
		return m;
	}
	uint64_t getNumberOfElements(){
		return data.size();
	}
};
}}
