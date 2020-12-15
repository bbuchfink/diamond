/*
The MIT License

Copyright (c) 2013-2019 Benjamin Buchfink

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef SPARSE_MATRIX_H_
#define SPARSE_MATRIX_H_

#include <vector>
#include <tuple>

struct SparseMatrix {

	typedef std::tuple<unsigned, unsigned, float> Triplet;

	SparseMatrix(std::vector<Triplet> &&v);
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