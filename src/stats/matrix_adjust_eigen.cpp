/****
DIAMOND protein aligner
Copyright (C) 2020-2021 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include "../lib/Eigen/Dense"
#include "../basic/value.h"

using namespace Eigen;

static const size_t N = TRUE_AA;
typedef double Float;
typedef Matrix<Float, N, N> MatrixN;
typedef Matrix<Float, 2 * N - 1, 2 * N - 1> Matrix2N;
typedef Matrix<Float, N, 1> VectorN;
typedef Matrix<Float, 2 * N - 1, 1> Vector2N;
typedef Matrix<Float, 2 * N, 1> Vector2Nx;
typedef Matrix<Float, 2, N * N> Vector2NN;
using Values = double[2];

static void ScaledSymmetricProductA(Matrix2N& W, const MatrixN& diagonal)
{
    int rowW, colW;   /* iteration indices over the rows and columns of W */
    int i, j;         /* iteration indices over characters in the alphabet */
    
    W.fill(0.0);
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            double dd;     /* an individual diagonal element */

            dd = diagonal(i, j);

            W(j, j) += dd;
            if (i > 0) {
                W(i + N - 1, j) += dd;
                W(i + N - 1, i + N - 1) += dd;
            }
        }
    }
}

static void MultiplyByA(Float beta, Vector2N& y, Float alpha, const MatrixN& x)
{
    int i, j;     /* iteration indices over characters in the alphabet */
    
    if (beta == 0.0)
        y.fill(0.0);
    else if (beta != 1.0) {
        /* rescale y */
        y *= beta;
    }

    VectorN sc = x.colwise().sum() * alpha, sr = x.rowwise().sum() * alpha;

    for (size_t i = 0; i < N; ++i)
        y(i) += sc(i);

    for (size_t i = 1; i < N; ++i)
        y(i + N - 1) += sr(i);
}

static void MultiplyByAtranspose(Float beta, MatrixN& y, Float alpha, const Vector2N& x)
{
    int i, j;     /* iteration indices over characters in the alphabet */
    int k;        /* index of a row of A transpose (a column of A); also
                      an index into y */

    if (beta == 0.0) {
        /* Initialize y to zero, without reading any elements from y */
        y.fill(0.0);
    }
    else if (beta != 1.0) {
        /* rescale y */
        y *= beta;
    }

    VectorN v = x.head<N>();
    v *= alpha;
    y.rowwise() += v.transpose();

    v = x.tail<N>();
    v[0] = 0.0;
    v *= alpha;
    y.colwise() += v;
}

static void ResidualsLinearConstraints(Vector2N& rA, const MatrixN& x, const VectorN& row_sums, const VectorN& col_sums)
{
    int i;             /* iteration index */

    rA.head<N>() = col_sums;
    rA.tail<N - 1>() = row_sums.tail<N - 1>();

    MultiplyByA(1.0, rA, -1.0, x);
}

static void DualResiduals(MatrixN& resids_x, const Vector2NN& grads, const Vector2Nx& z)
{        
    double eta = z[2 * N - 1];     /* dual variable for the relative
                                      entropy constraint */

    Matrix<Float, 1, N* N> v = -grads.row(0) + eta * grads.row(1);
    for (size_t i = 0; i < N; ++i)
        resids_x.col(i) = v.segment<N>(i * N);
    MultiplyByAtranspose(1.0, resids_x, 1.0, z.head<2*N-1>());
}

static void CalculateResiduals(double* rnorm,
    MatrixN& resids_x,
    const Vector2Nx& resids_z,
    const Values& values,
    double** grads,
    const double row_sums[],
    const double col_sums[],
    const double x[],
    const double z[],
    int constrain_rel_entropy,
    double relative_entropy)
{
    /* Euclidean norms of the primal and dual residuals */
    double norm_resids_z, norm_resids_x;

    DualResiduals(resids_x, alphsize, grads, z, constrain_rel_entropy);
    norm_resids_x = Nlm_EuclideanNorm(resids_x, alphsize * alphsize);

    ResidualsLinearConstraints(resids_z, alphsize, x, row_sums, col_sums);

    if (constrain_rel_entropy) {
        resids_z[2 * alphsize - 1] = relative_entropy - values[1];

        norm_resids_z = Nlm_EuclideanNorm(resids_z, 2 * alphsize);
    }
    else {
        norm_resids_z = Nlm_EuclideanNorm(resids_z, 2 * alphsize - 1);
    }
    *rnorm =
        sqrt(norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z);
}