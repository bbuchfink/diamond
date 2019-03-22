#include <algorithm>
#include "sparse_matrix.h"

using std::get;

SparseMatrix::SparseMatrix(std::vector<Triplet> &v)
{
	std::sort(v.begin(), v.end());
	const unsigned n = get<0>(v.back());
}