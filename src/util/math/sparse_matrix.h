#ifndef SPARSE_MATRIX_H_
#define SPARSE_MATRIX_H_

#include <iostream>
#include <vector>
#include <tuple>

struct SparseMatrix {

	typedef std::tuple<unsigned, unsigned, float> Triplet;

	SparseMatrix(std::vector<Triplet> &v);
	SparseMatrix(std::istream &is);
	SparseMatrix(unsigned rows);
	void print_stats();
	SparseMatrix transpose() const;
	friend SparseMatrix operator*(const SparseMatrix &A, const SparseMatrix &B);
	unsigned rows() const {
		return (unsigned)idx_.size();
	}

	unsigned nonzero;

private:

	static std::vector<Triplet> read_triplets(std::istream &is);
	static void multiply_worker(size_t i, size_t thread_id, const SparseMatrix *A, const SparseMatrix *B, SparseMatrix *C);

	std::vector<std::vector<unsigned>> idx_;
	std::vector<std::vector<float>> value_;

};

#endif