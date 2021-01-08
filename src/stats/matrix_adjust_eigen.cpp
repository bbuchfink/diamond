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
#include <type_traits>

using namespace Eigen;

static const size_t N = TRUE_AA;
typedef double Float;
typedef Matrix<Float, N, N> MatrixN;
typedef Matrix<Float, 2 * N, 2 * N> Matrix2Nx;
typedef Matrix<Float, N, 1> VectorN;
typedef Matrix<Float, 2 * N - 1, 1> Vector2N;
typedef Matrix<Float, 2 * N, 1> Vector2Nx;
typedef Matrix<Float, 1, N* N> VectorNN;
typedef Matrix<Float, 2, N * N> Vector2NN;
typedef decltype(Vector2Nx().head<2 * N - 1>()) Block2N;
using Values = double[2];

typedef struct ReNewtonSystem {
    Matrix2Nx W;               /**< A lower-triangular matrix
                                    representing a factorization of
                                    the (2,2) block, -J D^{-1} J^T, of
                                    the condensed linear system */
    MatrixN Dinv;             /**< The diagonal elements of the
                                    inverse of the necessarily
                                    diagonal (1,1) block of the linear
                                    system */
    VectorNN grad_re;          /**< the gradient of the
                                    relative-entropy constraint, if
                                    this constraint is used. */
} ReNewtonSystem;

static void ScaledSymmetricProductA(Matrix2Nx& W, const MatrixN& diagonal)
{    
    W.fill(0.0);    
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            const double dd = diagonal(i, j);

            W(j, j) += dd;
            if (i > 0) {
                W(i + N - 1, j) += dd;
                W(i + N - 1, i + N - 1) += dd;
            }
        }
    }
}

static void MultiplyByA(Float beta, Block2N y, Float alpha, const MatrixN& x)
{   
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

static void ResidualsLinearConstraints(Block2N rA, const MatrixN& x, const VectorN& row_sums, const VectorN& col_sums)
{
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
	Vector2Nx& resids_z,
	const Values& values,
	const Vector2NN& grads,
	const VectorN& row_sums,
	const VectorN& col_sums,
	const MatrixN& x,
	const Vector2Nx& z,
	double relative_entropy)
{
	/* Euclidean norms of the primal and dual residuals */
	double norm_resids_z, norm_resids_x;

	DualResiduals(resids_x, grads, z);
	norm_resids_x = resids_x.norm();

	ResidualsLinearConstraints(resids_z.head<2 * N - 1>(), x, row_sums, col_sums);
	resids_z[2 * N - 1] = relative_entropy - values[1];
	norm_resids_z = resids_z.norm();
	*rnorm = sqrt(norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z);
}

void Nlm_FactorLtriangPosDef(Matrix2Nx& A)
{
    int i, j, k;                /* iteration indices */
    double temp;                /* temporary variable for intermediate
                                   values in a computation */

    for (i = 0; i < 2*N; i++) {
        for (j = 0; j < i; j++) {
            temp = A(i, j);
            for (k = 0; k < j; k++) {
                temp -= A(i, k) * A(j, k);
            }
            A(i, j) = temp / A(j, j);
        }
        temp = A(i, i);
        for (k = 0; k < i; k++) {
            temp -= A(i, k) * A(i, k);
        }
        A(i,i) = sqrt(temp);
    }
}

static void FactorReNewtonSystem(ReNewtonSystem* newton_system,
    const MatrixN& x,
    const Vector2Nx& z,
    const Vector2NN& grads,
    MatrixN& workspace)
{
    int i;          /* iteration index */
    int n;          /* the length of x */
    int m;          /* the length of z */

    /* Pointers to fields in newton_systems; the names of the local
     * variables match the names of the fields. */
    Matrix2Nx& W = newton_system->W;
    MatrixN& Dinv = newton_system->Dinv;
    VectorNN& grad_re = newton_system->grad_re;

    n = N * N;
    m = 2 * N;

    /* The original system has the form
     *
     *     (D     J^T)
     *     (J     0  ).
     *
     * We block reduce the system to
     *
     *     (D    J^T          )
     *     (0    -J D^{-1} J^T).
     *
     * First we find the inverse of the diagonal matrix D. */

    double eta;             /* dual variable for the relative
                                   entropy constraint */
    eta = z[m - 1];

    Dinv = x;
    Dinv /= 1 - eta;

    /* Then we compute J D^{-1} J^T; First fill in the part that corresponds
     * to the linear constraints */
    ScaledSymmetricProductA(W, Dinv);

    /* Save the gradient of the relative entropy constraint. */
    grad_re = grads.row(1);

    /* Fill in the part of J D^{-1} J^T that corresponds to the relative
     * entropy constraint. */
    W(m - 1, m - 1) = 0.0;
    for (size_t i = 0; i < N; ++i) {
        workspace.col(i) = Dinv.col(i).cwiseProduct(grad_re.segment<N>(i * N).transpose());
        W(m - 1, m - 1) += workspace.col(i).dot(grad_re.segment<N>(i * N));
    }

    Vector2Nx r = W.row(m - 1);
    MultiplyByA(0.0, r.head<2 * N - 1>(), 1.0, workspace);

    /* Factor J D^{-1} J^T and save the result in W. */
    Nlm_FactorLtriangPosDef(W);
}
//
//static void SolveReNewtonSystem(const MatrixN& x, Vector2Nx& z,
//    const ReNewtonSystem* newton_system, MatrixN& workspace)
//{
//    int i;                     /* iteration index */
//    int n;                     /* the size of x */
//    int mA;                    /* the number of linear constraints */
//    int m;                     /* the size of z */
//
//    /* Local variables that represent fields of newton_system */
//    const Matrix2Nx& W = newton_system->W;
//    const MatrixN& Dinv = newton_system->Dinv;
//    const VectorNN& grad_re = newton_system->grad_re;
//    
//    n = N * N;
//    mA = 2 * N - 1;
//    m = mA + 1;
//
//    /* Apply the same block reduction to the right-hand side as was
//     * applied to the matrix:
//     *
//     *     rzhat = rz - J D^{-1} rx
//     */
//    workspace = x.cwiseProduct(Dinv);
//    MultiplyByA(1.0, z.head<2 * N - 1>(), -1.0, workspace);
//
//
//     for (i = 0; i < n; i++) {
//            z[m - 1] -= grad_re[i] * workspace[i];
//        }
//   
//
//    /* Solve for step in z, using the inverse of J D^{-1} J^T */
//    Nlm_SolveLtriangPosDef(z, m, W);
//
//    /* Backsolve for the step in x, using the newly-computed step in z.
//     *
//     *     x = D^{-1) (rx + J\T z)
//     */
//    if (constrain_rel_entropy) {
//        for (i = 0; i < n; i++) {
//            x[i] += grad_re[i] * z[m - 1];
//        }
//    }
//    MultiplyByAtranspose(1.0, x, alphsize, 1.0, z);
//
//    for (i = 0; i < n; i++) {
//        x[i] *= Dinv[i];
//    }
//}