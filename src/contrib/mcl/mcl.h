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
#include <limits>
#include <atomic>
#include "../util/system/system.h"
#include "../basic/config.h"
#include "../data/reference.h"
#include "../run/workflow.h"
#include "../util/io/consumer.h"
#include "../util/algo/algo.h"
#include "../basic/statistics.h"
#include "../util/log_stream.h"
#include "../dp/dp.h"
#include "../../cluster/cluster.h"
#include "sparse_matrix_stream.h"

namespace Workflow { namespace Cluster{
class MCL: public ClusteringAlgorithm {
private:
	using Weight = float;
	using Id = SparseMatrixStream<Weight>::Id;
	vector<Eigen::Triplet<float>> sparse_matrix_multiply(Eigen::SparseMatrix<float>* a, Eigen::SparseMatrix<float>* b , uint32_t iThr, uint32_t nThr);
	vector<Eigen::Triplet<float>> sparse_matrix_get_gamma(Eigen::SparseMatrix<float>* in, float r, uint32_t iThr, uint32_t nThr);
	float sparse_matrix_get_norm(Eigen::SparseMatrix<float>* in, uint32_t nThr);
	void print_stats(int64_t nElements, int64_t nComponents, int64_t nComponentsLt1, vector<int64_t>& sort_order, vector<vector<int64_t>>& indices, SparseMatrixStream<float>* ms);
	void get_exp(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r, uint32_t nThr);
	void get_exp(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r);
	void get_gamma(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r, uint32_t nThr);
	void get_gamma(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r);
	void markov_process(Eigen::SparseMatrix<float>* m, float inflation, float expansion, uint32_t max_iter, std::function<uint32_t()> getThreads);
	void markov_process(Eigen::MatrixXf* m, float inflation, float expansion, uint32_t max_iter);
	Eigen::SparseMatrix<float> get_sparse_matrix_and_clear(vector<int64_t>* order, vector<Eigen::Triplet<float>>* m, bool symmetric);
	Eigen::MatrixXf get_dense_matrix_and_clear(vector<int64_t>* order, vector<Eigen::Triplet<float>>* m, bool symmetric);
	std::atomic_ullong failed_to_converge = {0};
	std::atomic_ullong sparse_create_time = {0};
	std::atomic_ullong dense_create_time = {0};
	std::atomic_ullong sparse_exp_time = {0};
	std::atomic_ullong dense_int_exp_time = {0};
	std::atomic_ullong dense_gen_exp_time = {0};
	std::atomic_ullong sparse_gamma_time = {0};
	std::atomic_ullong dense_gamma_time = {0};
	std::atomic_ullong sparse_list_time = {0};
	std::atomic_ullong dense_list_time = {0};
public:
	~MCL(){};
	void run();
	string get_description();
	static string get_key(){
		return "mcl";
	}
};
}}
