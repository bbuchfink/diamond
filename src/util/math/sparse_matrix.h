#ifndef SPARSE_MATRIX_H_
#define SPARSE_MATRIX_H_

#include <vector>
#include <tuple>

struct SparseMatrix {

	typedef std::tuple<unsigned, unsigned, float> Triplet;

	SparseMatrix(std::vector<Triplet> &v);

private:

	std::vector<std::vector<unsigned>> idx_;
	std::vector<std::vector<float>> value_;

};

#endif