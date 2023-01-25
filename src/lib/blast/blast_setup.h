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
 * Author:  Tom Madden
 *
 */

/** @file blast_setup.h
 * Utilities initialize/setup BLAST.
 */

#ifndef __BLAST_SETUP__
#define __BLAST_SETUP__

#include "ncbi_std.h"
#include "blast_export.h"
#include "blast_def.h"
#include "blast_query_info.h"
#include "blast_options.h"

#include "blast_message.h"
#include "blast_stat.h"


#ifdef __cplusplus
extern "C" {
#endif

/** "Main" setup routine for BLAST. Calculates all information for BLAST search
 * that is dependent on the ASN.1 structures.
 * @todo FIXME: this function only filters query and sets up score block structure
 * @param program_number Type of BLAST program (0=blastn, ...). [in]
 * @param qsup_options options for query setup. [in]
 * @param scoring_options options for scoring. [in]
 * @param query_blk BLAST_SequenceBlk* for the query. [in]
 * @param query_info The query information block [in]
 * @param scale_factor Multiplier for cutoff and dropoff scores [in]
 * @param lookup_segments Start/stop locations for non-masked query 
 *                        segments [out]
 * @param mask masking locations. [out]
 * @param sbpp Contains scoring information. [out]
 * @param blast_message error or warning [out] 
 * @param get_path callback function to get matrix path [in]
 */
NCBI_XBLAST_EXPORT
Int2 BLAST_MainSetUp(EBlastProgramType program_number,
        const QuerySetUpOptions* qsup_options,
        const BlastScoringOptions* scoring_options,
        BLAST_SequenceBlk* query_blk,
        const BlastQueryInfo* query_info, 
        double scale_factor,
        BlastSeqLoc* *lookup_segments,
        BlastMaskLoc* *mask,
        BlastScoreBlk* *sbpp, 
        Blast_Message* *blast_message,
        GET_MATRIX_PATH get_path);

/** Blast_ScoreBlkKbpGappedCalc, fills the ScoreBlkPtr for a gapped search.  
 *      Should be moved to blast_stat.c in the future.
 * @param sbp Contains fields to be set, should not be NULL. [out]
 * @param scoring_options Scoring_options [in]
 * @param program Used to set fields on sbp [in]
 * @param query_info Query information containing context information [in]
 * @param error_return Pointer to structure for returning errors. [in][out]
 * @return Status.
 */
NCBI_XBLAST_EXPORT
Int2 Blast_ScoreBlkKbpGappedCalc(BlastScoreBlk * sbp,
                                 const BlastScoringOptions * scoring_options,
                                 EBlastProgramType program, 
                                 const BlastQueryInfo * query_info,
                                 Blast_Message** error_return);



NCBI_XBLAST_EXPORT
Int2 Blast_ScoreBlkMatrixInit(EBlastProgramType program_number, 
    const BlastScoringOptions* scoring_options,
    BlastScoreBlk* sbp,
    GET_MATRIX_PATH get_path);

/** Initializes the score block structure.
 * @param query_blk Query sequence(s) [in]
 * @param query_info Additional query information [in]
 * @param scoring_options Scoring options [in]
 * @param program_number BLAST program type [in]
 * @param sbpp Initialized score block [out]
 * @param scale_factor Matrix scaling factor for this search [in]
 * @param blast_message Error message [out]
 * @param get_path callback function to get matrix path [in]
 */
NCBI_XBLAST_EXPORT
Int2 BlastSetup_ScoreBlkInit(BLAST_SequenceBlk* query_blk, 
    const BlastQueryInfo* query_info, 
    const BlastScoringOptions* scoring_options, 
    EBlastProgramType program_number, 
    BlastScoreBlk* *sbpp, 
    double scale_factor, 
    Blast_Message* *blast_message,
    GET_MATRIX_PATH get_path);


/** Adjusts the mask locations coordinates to a sequence interval. Removes those
 * mask locations that do not intersect the interval. Can do this either for all 
 * queries or only for the first one.
 * @param mask Structure containing a mask location. [in] [out]
 * @param from Starting offset of a sequence interval [in]
 * @param to Ending offset of a sequence interval [in]
 */
NCBI_XBLAST_EXPORT
void
BlastSeqLoc_RestrictToInterval(BlastSeqLoc* *mask, Int4 from, Int4 to);



#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_SETUP__ */
