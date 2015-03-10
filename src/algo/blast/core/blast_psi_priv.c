#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] =
    "$Id: blast_psi_priv.c 366594 2012-06-15 15:08:33Z ucko $";
#endif /* SKIP_DOXYGEN_PROCESSING */
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
 * ===========================================================================
 *
 * Author:  Alejandro Schaffer, ported by Christiam Camacho
 *
 */

/** @file blast_psi_priv.c
 * Defintions for functions in private interface for Position Iterated BLAST 
 * API.
 */

#include "blast_psi_priv.h"
#include "ncbi_math.h"
#include "blast_util.h"

/****************************************************************************/
/* Use the following #define's to enable/disable functionality */

/* Taking gaps into account when constructing a PSSM was introduced in the 
 * 2001 paper "Improving the accuracy of PSI-BLAST protein database searches
 * with composition based-statistics and other refinements". This feature 
 * can be disabled by defining the PSI_IGNORE_GAPS_IN_COLUMNS symbol below */
/* #define PSI_IGNORE_GAPS_IN_COLUMNS */
/****************************************************************************/

/****************************************************************************/
/* Constants */
const double kPSINearIdentical = 0.94;
const double kPSIIdentical = 1.0;
const unsigned int kQueryIndex = 0;
const double kEpsilon = 0.0001;
const int kPSIScaleFactor = 200;
const double kPositScalingPercent = 0.05;
const Uint4 kPositScalingNumIterations = 10;

/****************************************************************************/

void**
_PSIAllocateMatrix(unsigned int ncols, unsigned int nrows, 
                   unsigned int data_type_sz)
{
    void** retval = NULL;
    unsigned int i = 0;

    retval = (void**) malloc(sizeof(void*) * ncols);
    if ( !retval ) {
        return NULL;
    }

    for (i = 0; i < ncols; i++) {
        retval[i] = (void*) calloc(nrows, data_type_sz);
        if ( !retval[i] ) {
            retval = _PSIDeallocateMatrix(retval, i);
            break;
        }
    }
    return retval;
}

void**
_PSIDeallocateMatrix(void** matrix, unsigned int ncols)
{
    unsigned int i = 0;

    if (!matrix)
        return NULL;

    for (i = 0; i < ncols; i++) {
        sfree(matrix[i]);
    }
    sfree(matrix);
    return NULL;
}

/** Implements the generic copy matrix functions. Prototypes must be defined
 * in the header file manually following the naming convention for 
 * _PSICopyMatrix_int
 */
#define DEFINE_COPY_MATRIX_FUNCTION(type)                           \
void _PSICopyMatrix_##type(type** dest, type** src,                 \
                          unsigned int ncols, unsigned int nrows)   \
{                                                                   \
    unsigned int i = 0;                                             \
    unsigned int j = 0;                                             \
                                                                    \
    ASSERT(dest);                                                   \
    ASSERT(src);                                                    \
                                                                    \
    for (i = 0; i < ncols; i++) {                                   \
        for (j = 0; j < nrows; j++) {                               \
            dest[i][j] = src[i][j];                                 \
        }                                                           \
    }                                                               \
}                                                                   \

DEFINE_COPY_MATRIX_FUNCTION(int)
DEFINE_COPY_MATRIX_FUNCTION(double)

/****************************************************************************/

/************** Validation routines *****************************************/

