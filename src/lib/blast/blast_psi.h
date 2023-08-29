/* $Id$
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

/** @file blast_psi.h
 * High level definitions and declarations for the PSSM engine of PSI-BLAST.
 */

#ifndef ALGO_BLAST_CORE___BLAST_PSI__H
#define ALGO_BLAST_CORE___BLAST_PSI__H

#include "ncbi_std.h"
#include "blast_export.h"
#include "blast_options.h"
#include "blast_stat.h"
#include <stdbool.h>


#ifdef __cplusplus
extern "C" {
#endif

/** Structure to describe the characteristics of a position in the multiple
 * sequence alignment data structure
 */
typedef struct PSIMsaCell {
    Uint1   letter;             /**< Preferred letter at this position, in
                                  ncbistdaa encoding */
    bool is_aligned;         /**< Is this letter part of the alignment? */
} PSIMsaCell;

/** Structure representing the dimensions of the multiple sequence alignment
 * data structure */
typedef struct PSIMsaDimensions {
    Uint4 query_length; /**< Length of the query */
    Uint4 num_seqs;     /**< Number of distinct sequences aligned with the 
                          query (does not include the query) */
} PSIMsaDimensions;

#ifdef DEBUG_PSSM_ENGINE
/** Sequence information structure */
typedef struct PSISeqInfo {
    int gi;                 /**< Sequence GI */
    double bit_score;       /**< Bit score of this sequence aligned with query*/
    double evalue;          /**< E-value of this sequence aligned with query*/
} PSISeqInfo;
#endif /* DEBUG_PSSM_ENGINE */

/** Multiple sequence alignment (msa) data structure containing the raw data 
 * needed by the PSSM engine to create a PSSM. By convention, the first row of
 * the data field contains the query sequence */
typedef struct PSIMsa {
    PSIMsaDimensions*   dimensions; /**< dimensions of the msa */
    PSIMsaCell**        data;       /**< actual data, dimensions are 
                                     (dimensions->num_seqs+1) by
                                     (dimensions->query_length) */
#ifdef DEBUG_PSSM_ENGINE
    PSISeqInfo*         seqinfo;    /** sequence information for a row of the
                                      multiple sequence alignment, length is
                                      dimensions->num_seqs+1 */
                                     
#endif /* DEBUG_PSSM_ENGINE */
} PSIMsa;

/** Allocates and initializes the multiple sequence alignment data structure
 * for use as input to the PSSM engine.
 * @param dimensions dimensions of multiple sequence alignment data structure
 * to allocate [in]
 * @return allocated PSIMsa structure or NULL if out of memory.
 */
NCBI_XBLAST_EXPORT
PSIMsa*
PSIMsaNew(const PSIMsaDimensions* dimensions);

/** Deallocates the PSIMsa structure
 * @param msa multiple sequence alignment structure to deallocate [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
PSIMsa*
PSIMsaFree(PSIMsa* msa);

/*******************************************************************
* Data structures for computing PSSM from using Conserved Domains
*
*/

/** Data needed for PSSM computation stored in MSA cell for single column in
 *  CD aligned to a position in the query */
typedef struct PSICdMsaCellData {

    double* wfreqs;                  /**< Frequencies for each residue in
                                          CD column */

    double iobsr;                    /**< Effective number of independent
                                          observations in a CD column */
} PSICdMsaCellData;

/** Alignment cell that represents one column of CD aligned to a position 
 *  in the query*/
typedef struct PSICdMsaCell {
    Uint1 is_aligned;          /**< Does this cell represent column aligned 
                                    to a CD */

    PSICdMsaCellData* data;    /**< Data needed for PSSM computation */
} PSICdMsaCell;


/** Data structure representing multiple alignemnt of CDs and query sequence
    along with data needed for PSSM computation */
typedef struct PSICdMsa {
    unsigned char* query;         /**< Query sequence as Ncbistdaa */
    PSIMsaDimensions* dimensions; /**< Query length and number of aligned cds */

    PSICdMsaCell **msa;          /**< Multiple alignment of CDs */
} PSICdMsa;


#ifdef DEBUG_PSSM_ENGINE
NCBI_XBLAST_EXPORT
void PrintMsa(const char* filename, const PSIMsa* msa);
NCBI_XBLAST_EXPORT
void PrintMsaFP(FILE* fp, const PSIMsa* msa);
#endif /* DEBUG_PSSM_ENGINE */

/** This is the main return value from the PSSM engine */
typedef struct PSIMatrix {
    Uint4   ncols;      /**< Number of columns in PSSM (query_length) */
    Uint4   nrows;      /**< Number of rows in PSSM (alphabet_size) */
    int**   pssm;       /**< Position-specific score matrix */
    double  lambda;     /**< Lambda Karlin-Altschul parameter */
    double  kappa;      /**< Kappa Karlin-Altschul parameter */
    double  h;          /**< H Karlin-Altschul parameter */
    double  ung_lambda; /**< Ungapped Lambda Karlin-Altschul parameter */
    double  ung_kappa;  /**< Ungapped Kappa Karlin-Altschul parameter */
    double  ung_h;      /**< Ungapped H Karlin-Altschul parameter */
} PSIMatrix;

/** Allocates a new PSIMatrix structure 
 * @param query_length number of columns allocated for the PSSM [in]
 * @param alphabet_size number of rows allocated for the PSSM [in]
 * @return pointer to allocated PSIMatrix structure or NULL if out of memory
 */
NCBI_XBLAST_EXPORT
PSIMatrix*
PSIMatrixNew(Uint4 query_length, Uint4 alphabet_size);

/** Deallocates the PSIMatrix structure passed in.
 * @param matrix structure to deallocate [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
PSIMatrix*
PSIMatrixFree(PSIMatrix* matrix);

/** Structure to allow requesting various diagnostics data to be collected by
 * PSSM engine */
typedef struct PSIDiagnosticsRequest {
    bool information_content;            /**< request information content */
    bool residue_frequencies;            /**< request observed residue
                                              frequencies */
    bool weighted_residue_frequencies;   /**< request observed weighted
                                              residue frequencies */
    bool frequency_ratios;               /**< request frequency ratios */
    bool gapless_column_weights;         /**< request gapless column weights
                                              */
    bool sigma;                          /**< request sigma */
    bool interval_sizes;                 /**< request interval sizes */
    bool num_matching_seqs;              /**< request number of matching 
                                              sequences */
    bool independent_observations;       /**< request number of independent
                                                 observations */

} PSIDiagnosticsRequest;

/** This structure contains the diagnostics information requested using the
 * PSIDiagnosticsRequest structure */
typedef struct PSIDiagnosticsResponse {
    double* information_content;           /**< position information content
                                             (query_length elements)*/
    Uint4** residue_freqs;                 /**< observed residue frequencies
                                             per position of the PSSM 
                                             (Dimensions are query_length by
                                             alphabet_size) */
    double** weighted_residue_freqs;       /**< Weighted observed residue
                                             frequencies per position of the
                                             PSSM. (Dimensions are query_length 
                                             by alphabet_size) */
    double** frequency_ratios;             /**< PSSM's frequency ratios
                                             (Dimensions are query_length by
                                             alphabet_size) */
    double* gapless_column_weights;        /**< Weights for columns without
                                             gaps (query_length elements) */
    double* sigma;                         /**< sigma (query_length elements) */
    Uint4* interval_sizes;                 /**< interval sizes of aligned
                                             regions (query_length elements) */
    Uint4* num_matching_seqs;              /**< number of matching sequences 
                                             per query position (query_length
                                             elements) */
    Uint4 query_length;                    /**< Specifies the number of
                                             positions in the PSSM */
    Uint4 alphabet_size;                   /**< Specifies length of alphabet */

    double* independent_observations;      /**< Effective number of
                                                observations per column */
} PSIDiagnosticsResponse;

/** Allocates a PSIDiagnosticsRequest structure, setting all fields to false
 * @return newly allocated structure or NULL in case of memory allocation
 * failure 
 */
NCBI_XBLAST_EXPORT
PSIDiagnosticsRequest* 
PSIDiagnosticsRequestNew(void);

/** Allocates a PSIDiagnosticsRequest structure, setting fields to their
 * default values for their use in the context of the PSI-BLAST application.
 * @param save_ascii_pssm corresponds to the command line argument to save the
 * PSSM in ASCII format [in]
 * @return newly allocated structure or NULL in case of memory allocation
 * failure 
 */
NCBI_XBLAST_EXPORT
PSIDiagnosticsRequest* 
PSIDiagnosticsRequestNewEx(bool save_ascii_pssm);

/** Deallocates the PSIDiagnosticsRequest structure passed in
 * @param diags_request structure to deallocate [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
PSIDiagnosticsRequest* 
PSIDiagnosticsRequestFree(PSIDiagnosticsRequest* diags_request);

/** Allocates a new PSI-BLAST diagnostics structure based on which fields of
 * the PSIDiagnosticsRequest structure are TRUE. Note: this is declared
 * here for consistency - this does not need to be called by client code of
 * this API, it is called in the PSICreatePssm* functions to allocate the 
 * diagnostics response structure.
 * @param query_length length of the query sequence [in]
 * @param alphabet_size length of the alphabet [in]
 * @param request diagnostics to retrieve from PSSM engine [in]
 * @return pointer to allocated PSIDiagnosticsResponse or NULL if dimensions or
 * request are NULL
 */
NCBI_XBLAST_EXPORT
PSIDiagnosticsResponse*
PSIDiagnosticsResponseNew(Uint4 query_length, Uint4 alphabet_size, 
                          const PSIDiagnosticsRequest* request);

/** Deallocates the PSIDiagnosticsResponse structure passed in.
 * @param diags structure to deallocate [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
PSIDiagnosticsResponse*
PSIDiagnosticsResponseFree(PSIDiagnosticsResponse* diags);

/****************************************************************************/

/** Main entry point to core PSSM engine to calculate the PSSM.
 * @param msap multiple sequence alignment data structure [in]
 * @param options options to the PSSM engine [in]
 * @param sbp BLAST score block structure [in|out]
 * @param pssm PSSM and statistical information (the latter is also returned 
 * in the sbp->kbp_gap_psi[0]) 
 * @return PSI_SUCCESS on success, otherwise one of the PSIERR_* constants
 */
NCBI_XBLAST_EXPORT
int
PSICreatePssm(const PSIMsa* msap,
              const PSIBlastOptions* options,
              BlastScoreBlk* sbp,
              PSIMatrix** pssm);

/** Main entry point to core PSSM engine which allows to request diagnostics
 * information.
 * @param msap multiple sequence alignment data structure [in]
 * @param options options to the PSSM engine [in]
 * @param sbp BLAST score block structure [in|out]
 * @param request diagnostics information request [in]
 * @param pssm PSSM and statistical information (the latter is also returned 
 * in the sbp->kbp_gap_psi[0]) [out]
 * @param diagnostics diagnostics information response, expects a pointer to an
 * uninitialized structure which will be populated with data requested in
 * requests [in|out]
 * @return PSI_SUCCESS on success, otherwise one of the PSIERR_* constants
 */
NCBI_XBLAST_EXPORT
int
PSICreatePssmWithDiagnostics(const PSIMsa* msap,
                             const PSIBlastOptions* options,
                             BlastScoreBlk* sbp,
                             const PSIDiagnosticsRequest* request,
                             PSIMatrix** pssm,
                             PSIDiagnosticsResponse** diagnostics);

/** Main entry point to core PSSM engine for computing CDD-based PSSMs
 * @param cd_msa information about CDs that match to query sequence [in]
 * @param options options to PSSM engine [in]
 * @param sbp BLAST score block structure [in|out]
 * @param request diagnostics information request [in]
 * @param pssm PSSM [out]
 * @param diagnostics diagnostics information response, expects a pointer to
 * uninitialized structure [in|out]
 */
int
PSICreatePssmFromCDD(const PSICdMsa* cd_msa,       
                     const PSIBlastOptions* options,
                     BlastScoreBlk* sbp, 
                     const PSIDiagnosticsRequest* request,
                     PSIMatrix** pssm,
                     PSIDiagnosticsResponse** diagnostics);



/** Top-level function to create a PSSM given a matrix of frequency ratios
 * and perform scaling on the resulting PSSM (i.e.: performs the last two 
 * stages of the algorithm)
 * Note that no diagnostics can be returned as those are calculated in earlier
 * stages of the algorithm.
 * @param query query sequence in ncbistdaa format, no sentinels needed [in]
 * @param query_length length of the query sequence [in]
 * @param sbp BLAST score block structure [in|out]
 * @param freq_ratios matrix of frequency ratios, dimensions are query_length
 * by BLASTAA_SIZE [in]
 * @param impala_scaling_factor scaling factor used in IMPALA-style scaling if
 * its value is NOT kPSSM_NoImpalaScaling (otherwise it performs standard
 * PSI-BLAST scaling) [in]
 * @param pssm PSSM and statistical information [in|out]
 * @return PSI_SUCCESS on success, otherwise one of the PSIERR_* constants
 * @todo FIXME change scalePosMatrix (blast_kappa.c) to use this function
 */
NCBI_XBLAST_EXPORT
int
PSICreatePssmFromFrequencyRatios(const Uint1* query,
                                 Uint4 query_length,
                                 BlastScoreBlk* sbp,
                                 double** freq_ratios,
                                 double impala_scaling_factor,
                                 PSIMatrix** pssm);

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__BLAST_PSI__H */

