/* ===========================================================================
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

/** @file nlm_linear_algebra.c
 * Basic matrix and vector operations
 *
 * @author E. Michael Gertz
 */

#include <math.h>
#include <stdlib.h>
//#include <algo/blast/core/ncbi_std.h>
#include "nlm_linear_algebra.h"


/* Documented in nlm_linear_algebra.h. */
double **
Nlm_DenseMatrixNew(size_t nrows,
                   size_t ncols)
{
    double ** mat;     /* the new matrix */

    mat = (double **) calloc(nrows, sizeof(double *));
    if (mat != NULL) {
        mat[0] = (double *) malloc((size_t) nrows *
                                   (size_t) ncols * sizeof(double));
        if (mat[0] != NULL) {
            for (size_t i = 1;  i < nrows;  i++) {
                mat[i] = &mat[0][i * ncols];
            }
        } else {
            free(mat);
            mat = NULL;
        }
    }
    return mat;
}


/* Documented in nlm_linear_algebra.h. */
double **
Nlm_LtriangMatrixNew(int n)
{
    int i;                      /* iteration index */
    double ** L;                /* the new, lower triangular matrix */
    size_t nelts;               /* the number of elements in
                                   the matrix */
    nelts = ((size_t) n * (n + 1))/2;

    L    = (double**) calloc(n, sizeof(double *));
    if (L != NULL) {
        L[0] = (double*) calloc(nelts, sizeof(double)); // was malloc
        if (L[0] != NULL) {
            for (i = 1;  i < n;  i++) {
                L[i] = L[i - 1] + i;
            }
        } else {
            free(L);
            L = NULL;
        }
    }
    return L;
}


/* Documented in nlm_linear_algebra.h. */
void
Nlm_DenseMatrixFree(double *** mat)
{
    if(*mat != NULL) {
        free((*mat)[0]);
        free(*mat);
    }
    *mat = NULL;
}


/* Documented in nlm_linear_algebra.h. */
int ** Nlm_Int4MatrixNew(int nrows, int ncols)
{
    int i;             /* iteration index */
    int ** mat;     /* the new matrix */

    mat = (int **) calloc(nrows, sizeof(int *));
    if (mat != NULL) {
        mat[0] = (int *) malloc((size_t) nrows *
                                   (size_t) ncols * sizeof(int));
        if (mat[0] != NULL) {
            for (i = 1;  i < nrows;  i++) {
                mat[i] = &mat[0][i * ncols];
            }
        } else {
            free(mat);
            mat = NULL;
        }
    }
    return mat;
}


/* Documented in nlm_linear_algebra.h. */
void
Nlm_Int4MatrixFree(int *** mat)
{
    if(*mat != NULL) {
        free((*mat)[0]);
        free(*mat);
    }
    *mat = NULL;
}


/* Documented in nlm_linear_algebra.h. */
void
Nlm_FactorLtriangPosDef(double ** A, int n)
{
    int i, j, k;                /* iteration indices */
    double temp;                /* temporary variable for intermediate
                                   values in a computation */

    for (i = 0;  i < n;  i++) {
        for (j = 0;  j < i;  j++) {
            temp = A[i][j];
            for (k = 0;  k < j;  k++) {
                temp -= A[i][k] * A[j][k];
            }
            A[i][j] = temp/A[j][j];
        }
        temp = A[i][i];
        for (k = 0;  k < i;  k++) {
            temp -= A[i][k] * A[i][k];
        }
        A[i][i] = sqrt(temp);
    }
}


/* Documented in nlm_linear_algebra.h. */
void Nlm_SolveLtriangPosDef(double x[], int n,
                            double ** L )
{
    int i, j;                   /* iteration indices */
    double temp;                /* temporary variable for intermediate
                                   values in a computation */

    /* At point x = b in the equation L L\T y = b */

    /* Forward solve; L z = b */
    for (i = 0;  i < n;  i++) {
        temp = x[i];
        for (j = 0;  j < i;  j++) {
            temp -= L[i][j] * x[j];
        }
        x[i] = temp/L[i][i];
    }
    /* Now x = z.  Back solve the system L\T y = z */
    for (j = n - 1;  j >= 0;  j--) {
        x[j] /= L[j][j];
        for (i = 0;  i < j;  i++) {
            x[i] -= L[j][i] * x[j];
        }
    }
    /* Now x = y, the solution to  L L\T y = b */
}


/* Documented in nlm_linear_algebra.h. */
double
Nlm_EuclideanNorm(const double v[], int n)
{
    double sum   = 1.0;   /* sum of squares of elements in v */
    double scale = 0.0;   /* a scale factor for the elements in v */
    int i;                /* iteration index */

    for (i = 0;  i < n;  i++) {
        if (v[i] != 0.0) {
            double absvi = fabs(v[i]);
            if (scale < absvi) {
                sum = 1.0 + sum * (scale/absvi) * (scale/absvi);
                scale = absvi;
            } else {
                sum += (absvi/scale) * (absvi/scale);
            }
        }
    }
    return scale * sqrt(sum);
}


/* Documented in nlm_linear_algebra.h. */
void Nlm_AddVectors(double y[], int n, double alpha, const double x[])
{
    int i;                     /* iteration index */

    for (i = 0; i < n; i++) {
        y[i] += alpha * x[i];
    }
}


/* Documented in nlm_linear_algebra.h. */
double
Nlm_StepBound(const double x[], int n, const double step_x[], double max)
{
    int i;                 /* iteration index */
    double alpha = max;    /* current largest permitted step */

    for (i = 0; i < n; i++) {
        double alpha_i;    /* a step to the boundary for the current i */

        alpha_i = -x[i] / step_x[i];
        if (alpha_i >= 0 && alpha_i < alpha) {
            alpha = alpha_i;
        }
    }
    return alpha;
}
