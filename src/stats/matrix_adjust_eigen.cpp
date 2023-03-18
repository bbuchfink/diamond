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

#include <type_traits>
#include "../lib/Eigen/Dense"
#include "../basic/value.h"
#include "../util/profiler.h"

// #define DYNAMIC

using namespace Eigen;
using std::endl;

//#define DEBUG_OUT(x) cout << (x) << endl
#define DEBUG_OUT(x)

static const size_t N = TRUE_AA;
typedef float Float;
const auto StorageOrder = RowMajor;
#ifdef DYNAMIC
typedef Matrix<Float, Dynamic, Dynamic, StorageOrder> MatrixN;
typedef Matrix<Float, Dynamic, Dynamic, StorageOrder> Matrix2Nx;
typedef Matrix<Float, Dynamic, 1> VectorN;
typedef Matrix<Float, Dynamic, 1> Vector2N;
typedef Matrix<Float, Dynamic, 1> Vector2Nx;
typedef Matrix<Float, 1, Dynamic> VectorNN;
typedef Matrix<Float, 2, Dynamic, StorageOrder> Vector2NN;
typedef decltype(Vector2Nx().head(Index())) Block2N;
#else
typedef Matrix<Float, N, N, StorageOrder> MatrixN;
typedef Matrix<Float, 2 * N, 2 * N, StorageOrder> Matrix2Nx;
typedef Matrix<Float, N, 1> VectorN;
typedef Matrix<Float, 2 * N - 1, 1> Vector2N;
typedef Matrix<Float, 2 * N, 1> Vector2Nx;
typedef Matrix<Float, 1, N* N> VectorNN;
typedef Matrix<Float, 2, N * N, StorageOrder> Vector2NN;
typedef decltype(Vector2Nx().head<2 * N - 1>()) Block2N;
#endif
using Values = Float[2];

typedef struct ReNewtonSystem {
#ifdef DYNAMIC
    ReNewtonSystem():
        W(2*N,2*N),
        Dinv(N,N),
        grad_re(N*N)
    {}
#endif
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
    Profiler prof("ScaledSymmetricProductA");
    W.fill(0.0);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            const Float dd = diagonal(i, j);

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
    Profiler prof("MultiplyByA");
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
    Profiler prof("MultiplyByAtranspose");
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
#ifdef DYNAMIC
    rA.head(N) = col_sums;
    rA.tail(N-1) = row_sums.tail(N-1);
#else
    rA.head<N>() = col_sums;
    rA.tail<N - 1>() = row_sums.tail<N - 1>();
#endif
    MultiplyByA(1.0, rA, -1.0, x);
    DEBUG_OUT("ResidualsLinearConstraints rA");
    DEBUG_OUT(rA);
}

static void DualResiduals(MatrixN& resids_x, const Vector2NN& grads, const Vector2Nx& z)
{        
    Float eta = z[2 * N - 1];     /* dual variable for the relative
                                      entropy constraint */

    VectorNN v = -grads.row(0) + eta * grads.row(1);
    for (size_t i = 0; i < N; ++i)
        resids_x.row(i) = v.segment<N>(i * N);
    MultiplyByAtranspose(1.0, resids_x, 1.0, z.head<2*N-1>());
    DEBUG_OUT("DualResiduals resids_x ");
    DEBUG_OUT(resids_x);
}

static void CalculateResiduals(Float* rnorm,
	MatrixN& resids_x,
	Vector2Nx& resids_z,
	const Values& values,
	const Vector2NN& grads,
	const VectorN& row_sums,
	const VectorN& col_sums,
	const MatrixN& x,
	const Vector2Nx& z,
	Float relative_entropy)
{
	/* Euclidean norms of the primal and dual residuals */
	Float norm_resids_z, norm_resids_x;

	DualResiduals(resids_x, grads, z);
	norm_resids_x = resids_x.norm();

#ifdef DYNAMIC
	ResidualsLinearConstraints(resids_z.head(2*N-1), x, row_sums, col_sums);
#else
    ResidualsLinearConstraints(resids_z.head<2 * N - 1>(), x, row_sums, col_sums);
#endif
	resids_z[2 * N - 1] = relative_entropy - values[1];
	norm_resids_z = resids_z.norm();
	*rnorm = sqrt(norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z);
    DEBUG_OUT("CalculateResiduals rnorm=");
    DEBUG_OUT(*rnorm);
}

static void FactorReNewtonSystem(ReNewtonSystem& newton_system,
    const MatrixN& x,
    const Vector2Nx& z,
    const Vector2NN& grads,
    MatrixN& workspace)
{
    Profiler prof("FactorReNewtonSystem");
    int n;          /* the length of x */
    int m;          /* the length of z */

    /* Pointers to fields in newton_systems; the names of the local
     * variables match the names of the fields. */
    Matrix2Nx& W = newton_system.W;
    MatrixN& Dinv = newton_system.Dinv;
    VectorNN& grad_re = newton_system.grad_re;

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

    Float eta;             /* dual variable for the relative
                                   entropy constraint */
    eta = z[m - 1];

    Dinv = x;
    Dinv /= 1 - eta;

    /* Then we compute J D^{-1} J^T; First fill in the part that corresponds
     * to the linear constraints */
    ScaledSymmetricProductA(W, Dinv);
    DEBUG_OUT("FactorReNewtonSystem W=");
    DEBUG_OUT(W);

    /* Save the gradient of the relative entropy constraint. */
    grad_re = grads.row(1);

    /* Fill in the part of J D^{-1} J^T that corresponds to the relative
     * entropy constraint. */
    W(m - 1, m - 1) = 0.0;
    for (size_t i = 0; i < N; ++i) {
#ifdef DYNAMIC
        workspace.row(i) = Dinv.row(i).cwiseProduct(grad_re.segment(i * N, N));
        W(m - 1, m - 1) += workspace.row(i).dot(grad_re.segment(i * N, N));
#else
        workspace.row(i) = Dinv.row(i).cwiseProduct(grad_re.segment<N>(i * N));
        W(m - 1, m - 1) += workspace.row(i).dot(grad_re.segment<N>(i * N));
#endif
    }

    Vector2Nx r = W.row(m - 1);
#ifdef DYNAMIC
    MultiplyByA(0.0, r.head(2*N-1), 1.0, workspace);
#else
    MultiplyByA(0.0, r.head<2 * N - 1>(), 1.0, workspace);
#endif
    W.row(m - 1) = r;
}

static void SolveReNewtonSystem(MatrixN& x, Vector2Nx& z, const ReNewtonSystem& newton_system, MatrixN& workspace) {
    Profiler prof("SolveReNewtonSystem");
    const Matrix2Nx& W = newton_system.W;
    const MatrixN& Dinv = newton_system.Dinv;
    const VectorNN& grad_re = newton_system.grad_re;

    const int m = 2 * N;

    /* Apply the same block reduction to the right-hand side as was
     * applied to the matrix:
     *
     *     rzhat = rz - J D^{-1} rx
     */
    workspace = x.cwiseProduct(Dinv);
#ifdef DYNAMIC
    MultiplyByA(1.0, z.head(2*N-1), -1.0, workspace);
#else
    MultiplyByA(1.0, z.head<2 * N - 1>(), -1.0, workspace);
#endif

    for (size_t i = 0; i < N; ++i)
        z[m - 1] -= grad_re.segment<N>(i * N).dot(workspace.row(i));

    /* Solve for step in z, using the inverse of J D^{-1} J^T */
    Profiler prof2("LLT");
    z = W.llt().solve(z);
    prof.finish();

    /* Backsolve for the step in x, using the newly-computed step in z.
     *
     *     x = D^{-1) (rx + J\T z)
     */
    for (size_t i = 0; i < N; ++i)
        x.row(i) += grad_re.segment<N>(i * N) * z[m - 1];

    MultiplyByAtranspose(1.0, x, 1.0, z.head<2 * N - 1>());
    x.array() *= Dinv.array();
}

static void EvaluateReFunctions(Values& values, Vector2NN& grads, const MatrixN& x, const MatrixN& q, const MatrixN& scores)
{
    Profiler prof("EvaluateReFunctions");
    values[0] = 0.0; values[1] = 0.0;

    MatrixN tmp = (x.array() / q.array()).log();
    for (size_t i = 0; i < N; ++i) {
        values[0] += x.row(i).dot(tmp.row(i));
        grads.row(0).segment<N>(i * N) = tmp.row(i);
        grads.row(0).segment<N>(i * N).array() += 1.0;
    }
    
    tmp.array() += scores.array();
    for (size_t i = 0; i < N; ++i) {
        values[1] += x.row(i).dot(tmp.row(i));
        grads.row(1).segment<N>(i * N) = tmp.row(i);
        grads.row(1).segment<N>(i * N).array() += 1.0;
    }
}

static void ComputeScoresFromProbs(MatrixN& scores,
    const MatrixN& target_freqs,
    const VectorN& row_freqs,
    const VectorN& col_freqs)
{
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            scores(i, j) = log(target_freqs(i, j) / (row_freqs[i] * col_freqs[j]));
        }
    }
}

static Float Nlm_StepBound(const MatrixN& x, const MatrixN& step_x, Float max)
{
    Float alpha = max;    /* current largest permitted step */

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; ++j) {
            Float alpha_i;    /* a step to the boundary for the current i */

            alpha_i = -x(i, j) / step_x(i, j);
            if (alpha_i >= 0 && alpha_i < alpha) {
                alpha = alpha_i;
            }
        }
    }
    return alpha;
}

static bool Blast_OptimizeTargetFrequencies(MatrixN& x,
    const MatrixN& q,
    const VectorN& row_sums,
    const VectorN& col_sums,
    Float relative_entropy,
    Float tol,
    int maxits)
{
    DEBUG_OUT("========================================================================================================");
#ifdef DYNAMIC
    Values values;   /* values of the nonlinear functions at this iterate */
    Vector2NN grads(2,N*N);     /* gradients of the nonlinear functions at this iterate */
    ReNewtonSystem newton_system;   /* factored matrix of the linear system to be solved at this iteration */
    Vector2Nx z(2*N);     /* dual variables (Lagrange multipliers) */
    z.fill(0.0);
    MatrixN resids_x(N,N);   /* dual residuals (gradient of Lagrangian) */
    Vector2Nx resids_z(2*N);   /* primal (constraint) residuals */
    Float rnorm;               /* norm of the residuals for the current iterate */
    MatrixN old_scores(N,N); /* a scoring matrix, with lambda = 1, generated from q, row_sums and col_sums */
    MatrixN workspace(N,N);  /* A vector for intermediate computations */
#else
    Values values;   /* values of the nonlinear functions at this iterate */
    Vector2NN grads;     /* gradients of the nonlinear functions at this iterate */
    ReNewtonSystem newton_system;   /* factored matrix of the linear system to be solved at this iteration */
    Vector2Nx z;     /* dual variables (Lagrange multipliers) */
    z.fill(0.0);
    MatrixN resids_x;   /* dual residuals (gradient of Lagrangian) */
    Vector2Nx resids_z;   /* primal (constraint) residuals */
    Float rnorm;               /* norm of the residuals for the current iterate */
    MatrixN old_scores; /* a scoring matrix, with lambda = 1, generated from q, row_sums and col_sums */
    MatrixN workspace;  /* A vector for intermediate computations */
#endif

    ComputeScoresFromProbs(old_scores, q, row_sums, col_sums);
    DEBUG_OUT(old_scores);

    /* Use q as the initial value for x */
    x = q;
    int its = 0;        /* Initialize the iteration count. Note that we may converge in zero iterations if the initial x is optimal. */
    while (its <= maxits) {
        /* Compute the residuals */
        EvaluateReFunctions(values, grads, x, q, old_scores);
        DEBUG_OUT("Values ");
        DEBUG_OUT(values[0]);
        DEBUG_OUT(values[1]);
        DEBUG_OUT("Grads");
        DEBUG_OUT(grads);

        CalculateResiduals(&rnorm, resids_x, resids_z, values,
            grads, row_sums, col_sums, x, z, relative_entropy);

        /* and check convergence; the test correctly handles the case
           in which rnorm is NaN (not a number). */
        if (!(rnorm > tol)) {
            /* We converged at the current iterate */
            break;
        }

        if (++its <= maxits) {
            FactorReNewtonSystem(newton_system, x, z, grads, workspace);
            SolveReNewtonSystem(resids_x, resids_z, newton_system, workspace);

            /* Calculate a value of alpha that ensure that x is positive */
            const Float alpha = Nlm_StepBound(x, resids_x, Float(1.0 / .95)) * Float(0.95);

            x += alpha * resids_x;
            z += alpha * resids_z;
        }

    }
    DEBUG_OUT(x);
    return its <= maxits && rnorm <= tol && z[2*N-1] < 1.0;
}

namespace Stats {

bool OptimizeTargetFrequencies(double* out, const double* joints_prob, const double* row_probs, const double* col_probs, double relative_entropy, double tol, int maxits) {
#ifdef DYNAMIC
    MatrixN x(N,N) , q(N,N);
    VectorN row_sums(N), col_sums(N);
#else
    MatrixN x, q;
    VectorN row_sums, col_sums;
#endif
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j)
            q(i, j) = (Float)joints_prob[i * N + j];
        row_sums[i] = (Float)row_probs[i];
        col_sums[i] = (Float)col_probs[i];
    }
    bool r = Blast_OptimizeTargetFrequencies(x, q, row_sums, col_sums, (Float)relative_entropy, (Float)tol, maxits);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            out[i * N + j] = x(i, j);
    return r;
}

}