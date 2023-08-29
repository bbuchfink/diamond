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

namespace Workflow { namespace Cluster{

void MCL::get_exp(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r){
	// TODO: at some r it may be more beneficial to diagnoalize in and only take the exponents of the eigenvalues
	high_resolution_clock::time_point t = high_resolution_clock::now();
	if(r - (int) r == 0){
		*out = *in;
		for(uint32_t i=1; i<r; i++){ 
			(*out) *= (*in);
		}
		//*out = in->pow(r)
		dense_int_exp_time += duration_cast<milliseconds>(high_resolution_clock::now() - t).count();
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
		dense_gen_exp_time += duration_cast<milliseconds>(high_resolution_clock::now() - t).count();
	}
}

void MCL::get_gamma(Eigen::SparseMatrix<float>* in, Eigen::SparseMatrix<float>* out, float r, uint32_t nThr) {
	high_resolution_clock::time_point t = high_resolution_clock::now();
	vector<Eigen::Triplet<float>> data;
	std::mutex m;
	auto mult = [&](const uint32_t iThr) {
		vector<Eigen::Triplet<float>> t_data = sparse_matrix_get_gamma(in, r, iThr, nThr);
		m.lock();
		data.insert(data.end(), t_data.begin(), t_data.end());
		m.unlock();
	};

	vector<thread> threads;
	for (uint32_t iThread = 0; iThread < nThr; iThread++) {
		threads.emplace_back(mult, iThread);
	}

	for (uint32_t iThread = 0; iThread < nThr; iThread++) {
		threads[iThread].join();
	}
	out->setZero();
	out->setFromTriplets(data.begin(), data.end(), [](const float&, const float &b) { return b; });
	*out = out->pruned(1.0, numeric_limits<float>::epsilon());
	sparse_gamma_time += duration_cast<milliseconds>(high_resolution_clock::now() - t).count();
}

void MCL::get_gamma(Eigen::MatrixXf* in, Eigen::MatrixXf* out, float r) {
	high_resolution_clock::time_point t = high_resolution_clock::now();
	// Note that Eigen matrices are column-major, so this is the most efficient way
	for (uint32_t icol = 0; icol < in->cols(); ++icol) {
		float colSum = 0.0f;
		for (uint32_t irow = 0; irow < in->rows(); ++irow) {
			colSum += pow(in->coeffRef(irow, icol), r);
		}
		for (uint32_t irow = 0; irow < in->rows(); ++irow) {
			out->coeffRef(irow, icol) = pow(in->coeffRef(irow, icol), r) / colSum;
		}
	}
	dense_gamma_time += duration_cast<milliseconds>(high_resolution_clock::now() - t).count();
}

void MCL::markov_process(Eigen::MatrixXf* m, float inflation, float expansion, uint32_t max_iter) {
	uint32_t iteration = 0;
	float diff_norm = numeric_limits<float>::max();
	Eigen::MatrixXf msquared(m->rows(), m->cols());
	Eigen::MatrixXf m_update(m->rows(), m->cols());
	get_gamma(m, m, 1); // This is to get a matrix of random walks on the graph -> TODO: find out if something else is more suitable
	while (iteration < max_iter && diff_norm > numeric_limits<float>::epsilon()) {
		get_exp(m, &msquared, expansion);
		get_gamma(&msquared, &m_update, inflation);
		*m -= m_update;
		diff_norm = m->norm();
		*m = m_update;
		iteration++;
	}
	if (iteration == max_iter) {
		failed_to_converge++;
	}
}

}}