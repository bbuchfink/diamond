/* $Id: nlm_linear_algebra.h 138123 2008-08-21 19:28:07Z camacho $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================*/

/**
 * @file nlm_linear_algebra.h
 * Declarations for several linear algebra routines
 *
 * @author E. Michael Gertz
 */

#ifndef __NLM_LINEAR_ALGEBRA__
#define __NLM_LINEAR_ALGEBRA__

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Create and return a new, dense matrix.  Elements of the matrix A
 * may be accessed as A[i][j]
 *
 * @param nrows     the number of rows for the new matrix.
 * @param ncols     the number of columns for the new matrix.
 */
double ** Nlm_DenseMatrixNew(size_t nrows, size_t ncols);


/**
 * Create and return a new, dense, lower-triangular matrix.  Elements
 * of the matrix A may be accessed as A[i][j] for i <= j.
 *
 * @param n         the dimension of the matrix.
 */
double ** Nlm_LtriangMatrixNew(int n);


/**
 * Free a matrix created by Nlm_DenseMatrixNew or
 * Nlm_LtriangMatrixNew.
 *
 * @param mat       the matrix to be freed
 * @return          always NULL
 */
void Nlm_DenseMatrixFree(double *** mat);


/**
 * Create and return a new Int4 matrix.  Elements of the matrix A
 * may be accessed as A[i][j]
 *
 * @param nrows     the number of rows for the new matrix.
 * @param ncols     the number of columns for the new matrix.
 */
int ** Nlm_Int4MatrixNew(int nrows, int ncols);


/**
 * Free a matrix created by Nlm_DenseMatrixNew or
 * Nlm_LtriangMatrixNew.
 *
 * @param mat       the matrix to be freed
 * @return          always NULL
 */
void Nlm_Int4MatrixFree(int *** mat);


/**
 * Accessing only the lower triangular elements of the symmetric,
 * positive definite matrix A, compute a lower triangular matrix L
 * such that A = L L^T (Cholesky factorization.)  Overwrite the lower
 * triangle of A with L.
 *
 * This routine may be used with the Nlm_SolveLtriangPosDef routine to
 * solve systems of equations.
 *
 * @param A         the lower triangle of a symmetric, positive-definite
 *                  matrix
 * @param n         the size of A
 */
void Nlm_FactorLtriangPosDef(double ** A, int n);


/**
 * Solve the linear system \f$ L L^T y = b \f$, where L is a non-singular
 * lower triangular matrix, usually computed using
 * the Nlm_FactorLtriangPosDef routine.
 *
 * @param x         on entry, the right hand size of the linear system
 *                  L L^T y = b; on exit the solution
 * @param n         the size of x
 * @param L         a non-singular lower triangular matrix
 */
void Nlm_SolveLtriangPosDef(double x[], int n, double ** L);


/**
 * Compute the Euclidean norm (2-norm) of a vector.
 *
 * This routine is based on the (freely available) BLAS routine dnrm2,
 * which handles the scale of the elements of v in a stable fashion.
 *
 * @param v      a vector
 * @param n      the length of v
 */
double Nlm_EuclideanNorm(const double v[], int n);


/**
 * Let y = y + alpha * x
 *
 * @param y         a vector
 * @param x         another vector
 * @param n         the length of x and y
 * @param alpha     a scale factor
 */
void Nlm_AddVectors(double y[], int n, double alpha,
                    const double x[]);


/**
 * Given a nonnegative vector x and a nonnegative scalar max, returns
 * the largest value in [0, max] for which x + alpha * step_x >= 0.
 *
 * @param x         a vector with nonnegative elements
 * @param step_x    another vector
 * @param n         the size of x and step_x
 * @param max       a nonnegative scalar
 */
double Nlm_StepBound(const double x[], int n,
                     const double step_x[], double max);

#ifdef __cplusplus
}
#endif

#endif
