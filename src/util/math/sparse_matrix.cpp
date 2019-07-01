#include <algorithm>
#include <string>
#include <map>
#include <assert.h>
#include <utility>
#include "../parallel/thread_pool.h"
#include "sparse_matrix.h"
#include "../../basic/config.h"
#include "../util.h"
#include "../data_structures/double_iterator.h"

using std::get;
using std::vector;
using std::string;
using std::map;
using std::cerr;
using std::endl;
using std::pair;

SparseMatrix::SparseMatrix(unsigned rows):
	idx_(rows),
	value_(rows)
{
}

SparseMatrix::SparseMatrix(std::vector<Triplet> &v)
{
	std::sort(v.begin(), v.end());
	const unsigned n = get<0>(v.back()) + 1;
	idx_.resize(n);
	value_.resize(n);
	for (const Triplet &t : v) {
		idx_[get<0>(t)].push_back(get<1>(t));
		value_[get<0>(t)].push_back(get<2>(t));
		++nonzero;
	}
}

std::vector<SparseMatrix::Triplet> SparseMatrix::read_triplets(std::istream &is) {
	map<string, unsigned> label2id;
	vector<Triplet> r;
	string n1, n2;
	float w;
	while (is >> n1 >> n2 >> w) {
		if (label2id.find(n1) == label2id.end())
			label2id[n1] = label2id.size();
		if (label2id.find(n2) == label2id.end())
			label2id[n2] = label2id.size();
		r.emplace_back(label2id[n1], label2id[n2], w);
	}
	return r;
}

SparseMatrix::SparseMatrix(std::istream &is):
	SparseMatrix(read_triplets(is))
{
}

void SparseMatrix::print_stats() {
	cerr << "Rows = " << rows() << endl;
	cerr << "Nonzero = " << nonzero << endl;
}

SparseMatrix SparseMatrix::transpose() const {
	SparseMatrix T(rows());
	for (unsigned i = 0; i < rows(); ++i)
		for (unsigned a = 0; a < idx_[i].size(); ++a) {
			unsigned j = idx_[i][a];
			T.idx_[j].push_back(i);
			T.value_[j].push_back(value_[i][a]);
		}

	typedef DoubleIterator<vector<unsigned>::iterator, vector<float>::iterator, unsigned, float> It;
	for (unsigned i = 0; i < T.rows(); ++i) {
		It begin{ T.idx_[i].begin(), T.value_[i].begin() }, end{ T.idx_[i].end(), T.value_[i].end() };
		//std::sort(begin, end);
	}
	return T;
}

void SparseMatrix::multiply_worker(size_t i, size_t thread_id, const SparseMatrix *A, const SparseMatrix *B, SparseMatrix *C) {
	for (unsigned j : A->idx_[i]) {

	}
}

SparseMatrix operator*(const SparseMatrix &A, const SparseMatrix &B) {
	assert(A.rows() == B.rows());
	SparseMatrix C(A.rows());
	Util::Parallel::scheduled_thread_pool_auto(config.threads_, A.rows(), SparseMatrix::multiply_worker, &A, &B, &C);
	return C;
}