#ifndef ALGO_BLAST_CORE___MATRIX_FREQ_RATIOS__H
#define ALGO_BLAST_CORE___MATRIX_FREQ_RATIOS__H

/*  $Id$
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
 * ===========================================================================
 *
 * Author:  Christiam Camacho
 *
 */

/** @file matrix_freq_ratios.h
 *  Interface to retrieve the frequency ratios for various scoring matrices.
 *
 *  See explanation in p 2996 of Nucleic Acids Research, 2001, Vol 29, No 14.
 */

#include "blast_encoding.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Stores the frequency ratios along with their bit scale factor */
typedef struct SFreqRatios {

    /** The actual frequency ratios */
    double**   data;

    /** Used to multiply the values in the above matrix to obtain scores in bit
     * units */
    int        bit_scale_factor;

} SFreqRatios;

/** Retrive the matrix's frequency ratios.
 * @param matrix_name Available options include:
 *          BLOSUM62
 *          BLOSUM62_20
 *          BLOSUM62_20A
 *          BLOSUM62_20B
 *          BLOSUM45
 *          BLOSUM80
 *          BLOSUM50
 *          BLOSUM90
 *          PAM30
 *          PAM70
 *          PAM250
 * @return NULL on error
 */
NCBI_XBLAST_EXPORT SFreqRatios*
_PSIMatrixFrequencyRatiosNew(const char* matrix_name);

/** Deallocate the frequency ratios structure */
NCBI_XBLAST_EXPORT SFreqRatios*
_PSIMatrixFrequencyRatiosFree(SFreqRatios* freq_ratios);

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__MATRIX_FREQ_RATIOS__H */
