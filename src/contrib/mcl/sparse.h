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

using std::runtime_error;

namespace Workflow { namespace Cluster{
		
vector<Eigen::Triplet<float>> MCL::sparse_matrix_multiply(Eigen::SparseMatrix<float>* a, Eigen::SparseMatrix<float>* b , uint32_t iThr, uint32_t nThr){
	int64_t n_cols = b->cols();
	int64_t n_rows = a->rows();
	std::vector<float> result_col(n_rows, 0.0);
	vector<Eigen::Triplet<float>> data;
	for (uint32_t j=0; j<n_cols; ++j) {
		if(j%nThr == iThr){
			fill(result_col.begin(), result_col.end(), 0.0f);
			for (Eigen::SparseMatrix<float>::InnerIterator  rhsIt(*b, j); rhsIt; ++rhsIt) {
				const float y = rhsIt.value();
				const int64_t k = rhsIt.row();
				for (Eigen::SparseMatrix<float>::InnerIterator lhsIt(*a, k); lhsIt; ++lhsIt) {
					const int64_t i = lhsIt.row();
					const float x = lhsIt.value();
					result_col[i] += x*y;
				}
			}
			for (int64_t i=0; i<n_rows; ++i) {
				if(abs(result_col[i]) > numeric_limits<float>::epsilon()) {
					data.emplace_back(i, j, result_col[i]);
				}
			}
		}
	}
	return data;
}

void MCL::get_exp(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r, uint32_t nThr){
	high_resolution_clock::time_point t = high_resolution_clock::now();
	if( r - (int) r == 0 ){
		// TODO: at some r it may be more beneficial to diagnoalize in and only take the exponents of the eigenvalues
		*out = *in;
		for(uint32_t i=1; i<r; i++){ 
			// Note: this is a workaround for a crash in Eigen
			vector<Eigen::Triplet<float>> data;
			std::mutex m;
			auto mult = [&](const uint32_t iThr){
				vector<Eigen::Triplet<float>> t_data = sparse_matrix_multiply(in, out, iThr, nThr);
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
		throw runtime_error("Eigen does not provide an eigenvalue solver for sparse matrices");
	}
	*out = out->pruned();
	sparse_exp_time += duration_cast<milliseconds>(high_resolution_clock::now() - t).count();
}

vector<Eigen::Triplet<float>> MCL::sparse_matrix_get_gamma(Eigen::SparseMatrix<float>* in, float r, uint32_t iThr, uint32_t nThr) {
	vector<Eigen::Triplet<float>> data;
	for (uint32_t k = 0; k < in->outerSize(); ++k) {
		if (k%nThr == iThr) {
			float colSum = 0.0f;
			for (Eigen::SparseMatrix<float>::InnerIterator it(*in, k); it; ++it) {
				colSum += pow(it.value(), r);
			}
			for (Eigen::SparseMatrix<float>::InnerIterator it(*in, k); it; ++it) {
				float val = pow(it.value(), r) / colSum;
				if (abs(val) > numeric_limits<float>::epsilon()) {
					data.emplace_back(it.row(), it.col(), val);
				}
			}
		}
	}
	return data;
}

float MCL::sparse_matrix_get_norm(Eigen::SparseMatrix<float>* in, uint32_t nThr) {
	vector<float> data(nThr);
	fill(data.begin(), data.end(), 0.0);
	auto norm = [&](const uint32_t iThr) {
		for (uint32_t k = 0; k < in->outerSize(); ++k) {
			if (k%nThr == iThr) {
				for (Eigen::SparseMatrix<float>::InnerIterator it(*in, k); it; ++it) {
					data[iThr] += pow(it.value(), 2.0f);
				}
			}
		}
	};

	vector<thread> threads;
	for (uint32_t iThread = 0; iThread < nThr; iThread++) {
		threads.emplace_back(norm, iThread);
	}

	for (uint32_t iThread = 0; iThread < nThr; iThread++) {
		threads[iThread].join();
	}
	return pow(accumulate(data.begin(), data.end(), 0.0f), 0.5f);
}

void MCL::markov_process(Eigen::SparseMatrix<float>* m, float inflation, float expansion, uint32_t max_iter, function<uint32_t()> getThreads) {
	uint32_t iteration = 0;
	float diff_norm = numeric_limits<float>::max();
	Eigen::SparseMatrix<float> msquared(m->rows(), m->cols());
	Eigen::SparseMatrix<float> m_update(m->rows(), m->cols());
	get_gamma(m, m, 1, getThreads()); // This is to get a matrix of random walks on the graph -> TODO: find out if something else is more suitable
	while (iteration < max_iter && diff_norm > numeric_limits<float>::epsilon()) {
		get_exp(m, &msquared, expansion, getThreads());
		get_gamma(&msquared, &m_update, inflation, getThreads());
		*m -= m_update;
		diff_norm = sparse_matrix_get_norm(m, getThreads());
		*m = m_update;
		iteration++;
	}
	if (iteration == max_iter) {
		failed_to_converge++;
	}
}

}}