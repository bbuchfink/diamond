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
#include "blast_psi.h"
#include "matrix_freq_ratios.h"
#include <stdbool.h>

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

/** Compact version of PSIMsa structure */
typedef struct _PSIPackedMsa {
    PSIMsaDimensions*   dimensions; /**< dimensions of the msa */
    _PSIPackedMsaCell** data;       /**< actual data, dimensions are 
                                     (dimensions->num_seqs+1) by
                                     (dimensions->query_length) */
    bool*            use_sequence;       /**< used to indicate whether a
                                              sequence should be used for
                                              further processing by the engine
                                              (length: num_seqs + 1) */
}_PSIPackedMsa;

/** Allocates and initializes the compact version of the PSIMsa structure
 * (makes a deep copy) for internal use by the PSSM engine.
 * @param msa multiple sequence alignment data structure provided by the user
 * [in]
 * @return newly allocated structure or NULL in case of memory allocation
 * failure
 */
NCBI_XBLAST_EXPORT 
_PSIPackedMsa*
_PSIPackedMsaNew(const PSIMsa* msa);

/** Deallocates the _PSIMsa data structure.
 * @param msa multiple sequence alignment data structure to deallocate [in] 
 * @return NULL
 */
NCBI_XBLAST_EXPORT 
_PSIPackedMsa*
_PSIPackedMsaFree(_PSIPackedMsa* msa);

/** Retrieve the number of aligned sequences in the compact multiple sequence
 * alignment
 * @param msa multiple sequence alignment data structure to deallocate [in] 
 */
NCBI_XBLAST_EXPORT 
unsigned int
_PSIPackedMsaGetNumberOfAlignedSeqs(const _PSIPackedMsa* msa);


/** Internal data structure to represent a position in the multiple sequence
 * alignment data structure @sa _PSIMsa */
typedef struct _PSIMsaCell {
    unsigned int letter:7;        /**< Preferred letter at this position */
    unsigned int is_aligned:1;    /**< Is this letter part of the alignment? */
    SSeqRange    extents;         /**< Extents of this aligned position */
} _PSIMsaCell;

/** Internal multiple alignment data structure used by the PSSM engine */
typedef struct _PSIMsa {
    PSIMsaDimensions*   dimensions;         /**< dimensions of field below */
    _PSIMsaCell**       cell;               /**< multiple sequence alignment
                                              matrix (dimensions: query_length 
                                              x num_seqs + 1) */
    Uint1*              query;              /**< query sequence (length:
                                              query_length) */
    Uint4**             residue_counts;     /**< matrix to keep track of the
                                              raw residue counts at each
                                              position of the multiple sequence
                                              alignment (dimensions: 
                                              query_length x alphabet_size) */
    Uint4               alphabet_size;      /**< number of elements in 
                                              alphabet */
    Uint4*              num_matching_seqs;  /**< number of sequences aligned at
                                              a given position in the multiple
                                              sequence alignment (length: 
                                              query_length). Corresponds to
                                              posit.h:posSearch->posCount */
} _PSIMsa;

/** Allocates and initializes the internal version of the PSIMsa structure
 * (makes a deep copy) for internal use by the PSSM engine.
 * @param packed_msa compact multiple sequence alignment data structure [in]
 * @param alphabet_size number of elements in the alphabet that makes up the
 * aligned characters in the multiple sequence alignment [in]
 * @return newly allocated structure or NULL in case of memory allocation
 * failure
 */
NCBI_XBLAST_EXPORT 
_PSIMsa*
_PSIMsaNew(const _PSIPackedMsa* packed_msa, Uint4 alphabet_size);

/** Deallocates the _PSIMsa data structure.
 * @param msa multiple sequence alignment data structure to deallocate [in] 
 * @return NULL
 */
NCBI_XBLAST_EXPORT 
_PSIMsa*
_PSIMsaFree(_PSIMsa* msa);

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

/** This structure keeps track of the regions aligned between the query
 * sequence and those that were not purged. It is used when calculating the
 * sequence weights (replaces posExtents in old code) */
typedef struct _PSIAlignedBlock {
    SSeqRange*  pos_extnt;      /**< Dynamically allocated array of size 
                                  query_length to keep track of the extents 
                                  of each aligned position */

    Uint4*      size;           /**< Dynamically allocated array of size 
                                  query_length that contains the size of the 
                                  intervals in the array above */
} _PSIAlignedBlock;

/** Allocates and initializes the _PSIAlignedBlock structure.
 * @param query_length length of the query sequence of the multiple
 * sequence alignment [in]
 * @return newly allocated structure or NULL in case of memory allocation
 * failure 
 */
NCBI_XBLAST_EXPORT 
_PSIAlignedBlock*
_PSIAlignedBlockNew(Uint4 query_length);

/** Deallocates the _PSIAlignedBlock structure.
 * @param aligned_blocks data structure to deallocate [in] 
 * @return NULL
 */
NCBI_XBLAST_EXPORT 
_PSIAlignedBlock*
_PSIAlignedBlockFree(_PSIAlignedBlock* aligned_blocks);

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

/** Allocates and initializes the _PSISequenceWeights structure.
 * @param dims structure containing the multiple sequence alignment dimensions
 * [in]
 * @param sbp score block structure initialized for the scoring system used
 * with the query sequence [in]
 * @return newly allocated structure or NULL in case of memory allocation
 * failure 
 */
NCBI_XBLAST_EXPORT 
_PSISequenceWeights*
_PSISequenceWeightsNew(const PSIMsaDimensions* dims, const BlastScoreBlk* sbp);

/** Deallocates the _PSISequenceWeights structure.
 * @param seq_weights data structure to deallocate [in] 
 * @return NULL
 */
NCBI_XBLAST_EXPORT 
_PSISequenceWeights*
_PSISequenceWeightsFree(_PSISequenceWeights* seq_weights);

/* Return values for internal PSI-BLAST functions */

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

/** Main function for keeping only those selected sequences for PSSM
 * construction (stage 2). After this function the multiple sequence alignment
 * data will not be modified.
 * @sa implementation of PSICreatePssmWithDiagnostics
 * @param msa multiple sequence alignment data structure [in]
 * @return PSIERR_BADPARAM if alignment is NULL; PSI_SUCCESS otherwise
 */
NCBI_XBLAST_EXPORT 
int 
_PSIPurgeBiasedSegments(_PSIPackedMsa* msa);

/** Main validation function for multiple sequence alignment structure. Should
 * be called after _PSIPurgeBiasedSegments.
 * @param msa multiple sequence alignment data structure [in]
 * @param ignored_unaligned_positions determines whether the unaligned
 * positions test should be performend or not [in]
 * @return One of the errors defined above if validation fails or bad
 * parameter is passed in, else PSI_SUCCESS
 */
NCBI_XBLAST_EXPORT 
int
_PSIValidateMSA(const _PSIMsa* msa, bool ignored_unaligned_positions);

/** Main function to compute aligned blocks' properties for each position 
 * within multiple alignment (stage 3) 
 * Corresponds to posit.c:posComputeExtents
 * @param msa multiple sequence alignment data structure [in]
 * @param aligned_block data structure describing the aligned blocks'
 * properties for each position of the multiple sequence alignment [out]
 * @return PSIERR_BADPARAM if arguments are NULL
 *         PSI_SUCCESS otherwise
 */
NCBI_XBLAST_EXPORT 
int
_PSIComputeAlignmentBlocks(const _PSIMsa* msa,
                           _PSIAlignedBlock* aligned_block);

/** Main function to calculate the sequence weights. Should be called with the
 * return value of PSIComputeAlignmentBlocks (stage 4)
 * Corresponds to posit.c:posComputeSequenceWeights
 * @param msa multiple sequence alignment data structure [in]
 * @param aligned_blocks data structure describing the aligned blocks'
 * properties for each position of the multiple sequence alignment [in]
 * @param nsg_compatibility_mode set to true to emulate the structure group's
 * use of PSSM engine in the cddumper application. By default should be FALSE
 * [in]
 * @param seq_weights data structure containing the data needed to compute the
 * sequence weights [out]
 * @return PSIERR_BADPARAM if arguments are NULL, PSIERR_OUTOFMEM in case of
 * memory allocation failure, PSIERR_BADSEQWEIGHTS if the sequence weights fail
 * to add up to 1.0, PSI_SUCCESS otherwise
 */
NCBI_XBLAST_EXPORT 
int
_PSIComputeSequenceWeights(const _PSIMsa* msa,
                           const _PSIAlignedBlock* aligned_blocks,
                           bool nsg_compatibility_mode,
                          _PSISequenceWeights* seq_weights);

/** Main function to calculate CD weights and combine weighted residue counts
 * from matched CDs
 * @param cd_msa multiple alignment of conserved domains data structure [in]
 * @param sbp BLAST score block [in]
 * @param options CDD-related options [in]
 * @param seq_weights data structure with CD frequencies [out]
 */
NCBI_XBLAST_EXPORT 
int
_PSIComputeFrequenciesFromCDs(const PSICdMsa* cd_msa,
                              BlastScoreBlk* sbp,
                              const PSIBlastOptions* options,
                              _PSISequenceWeights* seq_weights);
    

/** Main function to compute the PSSM's frequency ratios (stage 5).
 * Implements formula 2 in Nucleic Acids Research, 2001, Vol 29, No 14.
 * Corresponds to posit.c:posComputePseudoFreqs
 * @param msa multiple sequence alignment data structure [in]
 * @param seq_weights data structure containing the data needed to compute the
 * sequence weights [in]
 * @param sbp score block structure initialized for the scoring system used
 * with the query sequence [in]
 * @param aligned_blocks data structure describing the aligned blocks'
 * properties for each position of the multiple sequence alignment [in]
 * @param pseudo_count pseudo count constant [in]
 * @param nsg_compatibility_mode set to true to emulate the structure group's
 * use of PSSM engine in the cddumper application. By default should be FALSE
 * @param internal_pssm PSSM being computed [out]
 * @return PSIERR_BADPARAM if arguments are NULL, PSI_SUCCESS otherwise
 */
NCBI_XBLAST_EXPORT 
int
_PSIComputeFreqRatios(const _PSIMsa* msa,
                      const _PSISequenceWeights* seq_weights,
                      const BlastScoreBlk* sbp,
                      const _PSIAlignedBlock* aligned_blocks,
                      Int4 pseudo_count,
                      bool nsg_compatibility_mode,
                      _PSIInternalPssmData* internal_pssm);

/** Main function to compute CD-based PSSM's frequency ratios
 * @param cd_msa multiple alignment of CDs [in]
 * @param seq_weights contains weighted residue frequencies and effective number
 * of observations [in]
 * @param sbp initialized score block data structure [in]
 * @param pseudo_count pseudo count constant [in]
 * @param internal_pssm PSSM [out]
 * @return status
 */
NCBI_XBLAST_EXPORT 
int
_PSIComputeFreqRatiosFromCDs(const PSICdMsa* cd_msa,
                             const _PSISequenceWeights* seq_weights,
                             const BlastScoreBlk* sbp,
                             Int4 pseudo_count,
                             _PSIInternalPssmData* internal_pssm);


/** Converts the PSSM's frequency ratios obtained in the previous stage to a 
 * PSSM of scores. (stage 6) 
 * @param internal_pssm PSSM being computed [in|out]
 * @param query query sequence in ncbistdaa encoding. The length of this
 * sequence is read from internal_pssm->ncols [in]
 * @param sbp score block structure initialized for the scoring system used
 * with the query sequence [in]
 * @param std_probs array containing the standard residue probabilities [in]
 * @return PSIERR_BADPARAM if arguments are NULL, PSI_SUCCESS otherwise
 */
NCBI_XBLAST_EXPORT 
int
_PSIConvertFreqRatiosToPSSM(_PSIInternalPssmData* internal_pssm,
                            const Uint1* query,
                            const BlastScoreBlk* sbp,
                            const double* std_probs);

/** Scales the PSSM (stage 7)
 * @param query query sequence in ncbistdaa encoding. The length of this
 * sequence is read from internal_pssm->ncols [in]
 * @param std_probs array containing the standard background residue 
 * probabilities [in]
 * @param internal_pssm PSSM being computed [in|out]
 * @param sbp score block structure initialized for the scoring system used
 * with the query sequence [in|out]
 * @return PSIERR_BADPARAM if arguments are NULL, PSIERR_POSITIVEAVGSCORE if
 * the average score of the generated PSSM is positive, PSI_SUCCESS otherwise
 */
NCBI_XBLAST_EXPORT 
int
_PSIScaleMatrix(const Uint1* query,
                const double* std_probs,
                _PSIInternalPssmData* internal_pssm,
                BlastScoreBlk* sbp);

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

/** Marks the (start, stop] region corresponding to sequence seq_index in
 * alignment so that it is not further considered for PSSM calculation.
 * Note that the query sequence cannot be purged.
 * @param   msa multiple sequence alignment data  [in|out]
 * @param   seq_index index of the sequence of interested in alignment [in]
 * @param   start start of the region to remove [in]
 * @param   stop stop of the region to remove [in]
 * @return  PSIERR_BADPARAM if no alignment is given, or if seq_index or stop
 *          are invalid, PSI_SUCCESS otherwise
 */
NCBI_XBLAST_EXPORT 
int
_PSIPurgeAlignedRegion(_PSIPackedMsa* msa,
                       unsigned int seq_index,
                       unsigned int start,
                       unsigned int stop);

/** Counts the number of sequences matching the query per query position
 * (columns of the multiple alignment) as well as the number of residues
 * present in each position of the query.
 * Should be called after multiple alignment data has been purged from biased
 * sequences.
 * @param msa multiple sequence alignment structure [in|out]
 */
NCBI_XBLAST_EXPORT 
void
_PSIUpdatePositionCounts(_PSIMsa* msa);

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

/** Collects diagnostic information from the process of creating the PSSM 
 * @param msa multiple sequence alignment data structure [in]
 * @param aligned_block aligned regions' extents [in]
 * @param seq_weights sequence weights data structure [in]
 * @param internal_pssm structure containing PSSM's frequency ratios [in]
 * @param diagnostics output parameter [out]
 * @return PSI_SUCCESS on success, PSIERR_OUTOFMEM if memory allocation fails
 * or PSIERR_BADPARAM if any of its arguments is NULL
 */
NCBI_XBLAST_EXPORT 
int
_PSISaveDiagnostics(const _PSIMsa* msa,
                    const _PSIAlignedBlock* aligned_block,
                    const _PSISequenceWeights* seq_weights,
                    const _PSIInternalPssmData* internal_pssm,
                    PSIDiagnosticsResponse* diagnostics);

/** Collects diagnostic information from the process of creating the CDD-based
 * PSSM 
 * @param cd_msa multiple alignment of CDs data structure [in]
 * @param seq_weights sequence weights data structure [in]
 * @param internal_pssm structure containing PSSM's frequency ratios [in]
 * @param diagnostics output parameter [out]
 * @return PSI_SUCCESS on success, PSIERR_OUTOFMEM if memory allocation fails
 * or PSIERR_BADPARAM if any of its arguments is NULL
 */
int
_PSISaveCDDiagnostics(const PSICdMsa* msa,
                      const _PSISequenceWeights* seq_weights,
                      const _PSIInternalPssmData* internal_pssm,
                      PSIDiagnosticsResponse* diagnostics);


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

/** Enable NCBI structure group customization to discard the query sequence,
 * as this really isn't the result of a PSI-BLAST iteration, but rather an
 * artificial consensus sequence of the multiple sequence alignment
 * constructed by them. This should be called after _PSIPurgeBiasedSegments.
 */
NCBI_XBLAST_EXPORT 
void
_PSIStructureGroupCustomization(_PSIMsa* msa);

/** Structure group validation function for multiple sequence alignment 
 * structure. Should be called after _PSIStructureGroupCustomization.
 * @param msa multiple sequence alignment data structure [in]
 * @return One of the errors defined above if validation fails or bad
 * parameter is passed in, else PSI_SUCCESS
 */
NCBI_XBLAST_EXPORT 
int
_PSIValidateMSA_StructureGroup(const _PSIMsa* msa);

/** Validation of multiple alignment of conserved domains structure
 * @param cd_msa multiple alignment of CDs [in]
 * @param alphabet_size alphabet size [in]
 * @return One of the errors defined above if validation fails or bad
 * parameter is passed in, else PSI_SUCCESS
 */
NCBI_XBLAST_EXPORT 
int
_PSIValidateCdMSA(const PSICdMsa* cd_msa, Uint4 alphabet_size);

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__BLAST_PSI_PRIV__H */
