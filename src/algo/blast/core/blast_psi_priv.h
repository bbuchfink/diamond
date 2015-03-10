/*  $Id: blast_psi_priv.h 347365 2011-12-16 13:17:19Z ivanov $
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
 * Author:  Alejandro Schaffer, ported by Christiam Camacho
 *
 */

/** @file blast_psi_priv.h
 * Private interface for Position Iterated BLAST API, contains the
 * PSSM generation engine.
 *
 * <pre>
 * Calculating PSSMs from Seq-aligns is a multi-stage process. These stages
 * include:
 * 1) Processing the Seq-align
 *      Examine alignment and extract information about aligned characters,
 *      performed at the API level
 * 2) Purge biased sequences: construct M multiple sequence alignment as
 *      described in page 3395[1] - performed at the core level; custom
 *      selection of sequences should be performed at the API level.
 * 3) Compute extents of the alignment: M sub C as described in page 3395[1]
 * 4) Compute sequence weights
 * 5) Compute residue frequencies
 * 6) Convert residue frequencies to PSSM
 * 7) Scale the resulting PSSM
 * </pre>
 */

#ifndef ALGO_BLAST_CORE___BLAST_PSI_PRIV__H
#define ALGO_BLAST_CORE___BLAST_PSI_PRIV__H

#include "ncbi_std.h"
#include "blast_stat.h"

#ifdef __cplusplus
extern "C" {
#endif

/****************************************************************************/
/* Extern declarations for constants (defined in blast_psi_priv.c) */

/** Percent identity threshold for discarding near-identical matches */
NCBI_XBLAST_EXPORT 
extern const double kPSINearIdentical;

/** Percent identity threshold for discarding identical matches */
NCBI_XBLAST_EXPORT 
extern const double kPSIIdentical;

/** Index into multiple sequence alignment structure for the query sequence */
NCBI_XBLAST_EXPORT 
extern const unsigned int kQueryIndex;

/** Small constant to test against 0 */
NCBI_XBLAST_EXPORT 
extern const double kEpsilon;

/** Successor to POSIT_SCALE_FACTOR  */
NCBI_XBLAST_EXPORT 
extern const int kPSIScaleFactor;

/** Constant used in scaling PSSM routines: Successor to POSIT_PERCENT */
NCBI_XBLAST_EXPORT 
extern const double kPositScalingPercent;
/** Constant used in scaling PSSM routines: Successor to POSIT_NUM_ITERATIONS */
NCBI_XBLAST_EXPORT 
extern const Uint4 kPositScalingNumIterations;

/****************************************************************************/
/* Matrix utility functions */

/** Generic 2 dimensional matrix allocator.
 * Allocates a ncols by nrows matrix with cells of size data_type_sz. Must be
 * freed using x_DeallocateMatrix
 * @param   ncols number of columns in matrix [in]
 * @param   nrows number of rows in matrix [in]
 * @param   data_type_sz size of the data type (in bytes) to allocate for each
 *          element in the matrix [in]
 * @return pointer to allocated memory or NULL in case of failure
 */
NCBI_XBLAST_EXPORT 
void**
_PSIAllocateMatrix(unsigned int ncols, unsigned int nrows, 
                   unsigned int data_type_sz);

/** Generic 2 dimensional matrix deallocator.
 * Deallocates the memory allocated by x_AllocateMatrix
 * @param matrix matrix to deallocate   [in]
 * @param ncols number of columns in the matrix [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT 
void**
_PSIDeallocateMatrix(void** matrix, unsigned int ncols);

/** Copies src matrix into dest matrix, both of which must be int matrices with
 * dimensions ncols by nrows
 * @param dest Destination matrix           [out]
 * @param src Source matrix                 [in]
 * @param ncols Number of columns to copy   [in]
 * @param nrows Number of rows to copy      [in]
 */
NCBI_XBLAST_EXPORT 
void
_PSICopyMatrix_int(int** dest, int** src,
                   unsigned int ncols, unsigned int nrows);

/** Copies src matrix into dest matrix, both of which must be double matrices 
 * with dimensions ncols by nrows
 * @param dest Destination matrix           [out]
 * @param src Source matrix                 [in]
 * @param ncols Number of columns to copy   [in]
 * @param nrows Number of rows to copy      [in]
 */
NCBI_XBLAST_EXPORT 
void
_PSICopyMatrix_double(double** dest, double** src,
                      unsigned int ncols, unsigned int nrows);

/****************************************************************************/
/* Structure declarations */

/** Compact version of the PSIMsaCell structure */
typedef struct _PSIPackedMsaCell {
    unsigned int letter:7;      /**< Preferred letter at this position, in
                                  ncbistdaa encoding */
    unsigned int is_aligned:1;  /**< Is this letter part of the alignment? */
} _PSIPackedMsaCell;

/** Internal representation of a PSSM in various stages of its creation and its
 * dimensions. */
typedef struct _PSIInternalPssmData {
    Uint4       ncols;          /**< number of columns (query_length) */
    Uint4       nrows;          /**< number of rows (alphabet_size) */
    int**       pssm;           /**< PSSM (scores) */
    int**       scaled_pssm;    /**< scaled PSSM (scores) */
    double**    freq_ratios;    /**< frequency ratios */
    double*     pseudocounts;   /**< pseudocount constant for each column */
} _PSIInternalPssmData;

/** Allocates a new _PSIInternalPssmData structure.
 * @param query_length number of columns for the PSSM [in]
 * @param alphabet_size number of rows for the PSSM [in]
 * @return newly allocated structure or NULL in case of memory allocation
 * failure 
 */
NCBI_XBLAST_EXPORT 
_PSIInternalPssmData*
_PSIInternalPssmDataNew(Uint4 query_length, Uint4 alphabet_size);

/** Deallocates the _PSIInternalPssmData structure.
 * @param pssm data structure to deallocate [in] 
 * @return NULL
 */
NCBI_XBLAST_EXPORT 
_PSIInternalPssmData*
_PSIInternalPssmDataFree(_PSIInternalPssmData* pssm);

/** Internal data structure to keep computed sequence weights */
typedef struct _PSISequenceWeights {
    /* Some notes on how elements of this structures corresponds to elements 
       of structures defined in posit.h:

    _PSISequenceWeights->match_weights is the same as posSearch->posMatchWeights
    _PSISequenceWeights->gapless_column_weights is the same as 
          posSearch->posGaplessColumnWeights
    _PSISequenceWeights->norm_seq_weights is the same as posSearch->posA
    _PSISequenceWeights->row_sigma is the same as posSearch->posRowSigma
    */
    double** match_weights;     /**< weighted observed residue frequencies (f_i
                                  in 2001 paper). (dimensions: query_length 
                                  x BlastScoreBlk's alphabet_size) */
    Uint4 match_weights_size;   /**< kept for help deallocate the field above */

    double* norm_seq_weights;   /**< Stores the normalized sequence weights
                                  (length: num_seqs + 1) */
    double* row_sigma;  /**< array of length num_seqs + 1 */
    /* Sigma: number of different characters occurring in matches within a
     * multi-alignment block - FIXME: why is it a double? */
    double* sigma;      /**< array of length query_length */

    double* std_prob;   /**< standard amino acid probabilities */

    /* This fields is required for important diagnostic output, they are
     * copied into diagnostics structure */
    double* gapless_column_weights; /**< FIXME */
    int** posDistinctDistrib; /**< For position i, how many positions in its block
                               have j distinct letters.  Copied from
                               posit.h:posSearchItems.posDistinctDistrib */
    Uint4 posDistinctDistrib_size; /**< Kept to deallocate field above. */
    int* posNumParticipating; /**< number of sequences at each position.  Copied from
                                posit.h:posSearchItems.posNumParticipating  */
    double* independent_observations; /**< Number of independent sequences
                                       per column */

} _PSISequenceWeights;

/** Successful operation */
#define PSI_SUCCESS             (0)
/** Bad parameter used in function */
#define PSIERR_BADPARAM         (-1)
/** Out of memory */
#define PSIERR_OUTOFMEM         (-2)   
/** Sequence weights do not add to 1 */
#define PSIERR_BADSEQWEIGHTS    (-3)   
/** No frequency ratios were found for the given scoring matrix */
#define PSIERR_NOFREQRATIOS     (-4)   
/** Positive average score found when scaling matrix */
#define PSIERR_POSITIVEAVGSCORE (-5)   
/** After purge stage of PSSM creation, no sequences are left */
#define PSIERR_NOALIGNEDSEQS    (-6)
/** GAP residue found in query sequence */
#define PSIERR_GAPINQUERY       (-7)
/** Found an entire column with no participating sequences */
#define PSIERR_UNALIGNEDCOLUMN  (-8)
/** Found an entire column full of GAP residues */
#define PSIERR_COLUMNOFGAPS     (-9)
/** Found flanking gap at start of alignment */
#define PSIERR_STARTINGGAP      (-10)
/** Found flanking gap at end of alignment */
#define PSIERR_ENDINGGAP        (-11)
/** Errors in conserved domain profile */
#define PSIERR_BADPROFILE       (-12)
/** Unknown error */
#define PSIERR_UNKNOWN          (-255)

/****************************************************************************/
/* Function prototypes for the various stages of the PSSM generation engine */

/****************************************************************************/
/* Function prototypes for auxiliary functions for the stages above */


/** Updates the Karlin-Altschul parameters based on the query sequence and 
 * PSSM's score frequencies. Port of blastool.c's updateLambdaK
 * @param pssm PSSM [in]
 * @param query query sequence in ncbistdaa encoding [in]
 * @param query_length length of the query sequence above [in]
 * @param std_probs array containing the standard background residue 
 * probabilities [in]
 * @param sbp Score block structure where the calculated lambda and K will be
 * returned [in|out]
 */
NCBI_XBLAST_EXPORT 
void
_PSIUpdateLambdaK(const int** pssm,
                  const Uint1* query,
                  Uint4 query_length,
                  const double* std_probs,
                  BlastScoreBlk* sbp);

/** Provides a similar function to _PSIScaleMatrix but it performs the scaling
 * as IMPALA did, i.e.: allowing the specification of a scaling factor and when
 * calculating the score probabilities, the query length includes 'X' residues.
 * @todo Ideally all scaling code should be refactored so that it is
 * consolidated, eliminating the need for blast_posit.[hc]. Please note that
 * blast_kappa.c's scalePosMatrix also does something very similar.
 * @todo remove std_probs as it's not used
 */
NCBI_XBLAST_EXPORT 
int
_IMPALAScaleMatrix(const Uint1* query, const double* std_probs,
                   _PSIInternalPssmData* internal_pssm, 
                   BlastScoreBlk* sbp,
                   double scaling_factor);
/** Calculates the length of the sequence without including any 'X' residues.
 * used in kappa.c
 * @param seq sequence to examine [in]
 * @param length length of the sequence above [in]
 * @return number of non-X residues in the sequence
 */
NCBI_XBLAST_EXPORT 
Uint4
_PSISequenceLengthWithoutX(const Uint1* seq, Uint4 length);

/** Compute the probabilities for each score in the PSSM.
 * This is only valid for protein sequences.
 * FIXME: Should this be moved to blast_stat.[hc]?
 * used in kappa.c in notposfillSfp()
 * @param pssm PSSM for which to compute the score probabilities [in]
 * @param query query sequence for the PSSM above in ncbistdaa encoding [in]
 * @param query_length length of the query sequence above [in]
 * @param std_probs array containing the standard background residue 
 * probabilities [in]
 * @param sbp score block structure initialized for the scoring system used
 * with the query sequence [in]
 * @return structure containing the score frequencies, or NULL in case of error
 */
NCBI_XBLAST_EXPORT 
Blast_ScoreFreq*
_PSIComputeScoreProbabilities(const int** pssm,
                              const Uint1* query,
                              Uint4 query_length,
                              const double* std_probs,
                              const BlastScoreBlk* sbp);

/** Calculates the information content from the scoring matrix
 * @param score_mat alphabet by alphabet_sz matrix of scores (const) [in]
 * @param std_prob standard residue probabilities [in]
 * @param query query sequence [in]
 * @param query_length length of the query [in]
 * @param alphabet_sz length of the alphabet used by the query [in]
 * @param lambda lambda parameter [in] FIXME documentation
 * @return array of length query_length containing the information content per
 * query position or NULL on error (e.g.: out-of-memory or NULL parameters)
 */
NCBI_XBLAST_EXPORT 
double*
_PSICalculateInformationContentFromScoreMatrix(
    Int4** score_mat,
    const double* std_prob,
    const Uint1* query,
    Uint4 query_length,
    Uint4 alphabet_sz,
    double lambda);

/** Calculates the information content from the residue frequencies calculated
 * in stage 5 of the PSSM creation algorithm 
 * Corresponds to posit.c:posFreqsToInformation
 * @sa _PSIComputeFreqRatios: stage 5
 * @param freq_ratios matrix of frequency ratios (dimensions: query_length x 
 * alphabet_sz) (const) [in]
 * @param std_prob standard residue probabilities [in]
 * @param query_length length of the query [in]
 * @param alphabet_sz length of the alphabet used by the query [in]
 * @return array of length query_length containing the information content per
 * query position or NULL on error (e.g.: out-of-memory or NULL parameters)
 */
NCBI_XBLAST_EXPORT 
double*
_PSICalculateInformationContentFromFreqRatios(
    double** freq_ratios,
    const double* std_prob,
    Uint4 query_length,
    Uint4 alphabet_sz);

#ifdef DEBUG_PSSM_ENGINE
void __printMsa(const char* filename, const _PSIPackedMsa* msa);
void __printMsaFP(FILE* fp, const _PSIPackedMsa* msa);
#endif /* DEBUG_PSSM_ENGINE */

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__BLAST_PSI_PRIV__H */
