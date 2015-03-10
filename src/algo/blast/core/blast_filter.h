/* $Id: blast_filter.h 191335 2010-05-12 12:50:24Z madden $
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
 * Author:  Ilya Dondoshansky
 *
 */

/** @file blast_filter.h
 * BLAST filtering functions. @todo FIXME: contains more than filtering 
 * functions, combine with blast_dust.h?
 */

#ifndef ALGO_BLAST_CORE__BLAST_FILTER__H
#define ALGO_BLAST_CORE__BLAST_FILTER__H

#include "ncbi_std.h"
#include "blast_def.h"
#include "blast_program.h"
#include "blast_query_info.h"
#include "blast_message.h"
#include "blast_options.h"

#ifdef __cplusplus
extern "C" {
#endif

/** BLASTNA element used to mask bases in BLAST */
NCBI_XBLAST_EXPORT
extern const Uint1 kNuclMask;
/** NCBISTDAA element used to mask residues in BLAST */
NCBI_XBLAST_EXPORT
extern const Uint1 kProtMask;

/** Repeats filtering default options. */
#define REPEATS_SEARCH_EVALUE 0.1       /**< Default e-value threshold, keep for C toolkit */
#define REPEATS_SEARCH_MINSCORE 26       /**< Default score cutoff */
#define REPEATS_SEARCH_PENALTY -1       /**< Default mismatch penalty */
#define REPEATS_SEARCH_REWARD 1       /**< Default match reward */
#define REPEATS_SEARCH_GAP_OPEN 2       /**< Default gap opening cost */
#define REPEATS_SEARCH_GAP_EXTEND 1     /**< Default gap extension cost */
#define REPEATS_SEARCH_WORD_SIZE 11     /**< Default word size */
#define REPEATS_SEARCH_XDROP_UNGAPPED 40/**< Default X-dropoff for ungapped 
                                           extension */
#define REPEATS_SEARCH_XDROP_FINAL 90   /**< Default X-dropoff for gapped 
                                           extension with traceback */
#define REPEATS_SEARCH_FILTER_STRING "F"/**< Default filter string - 
                                           no filtering */

/** Largest gap allowed to be filled between repeat mask intervals */
#define REPEAT_MASK_LINK_VALUE 5

/** Create and initialize a new sequence interval.
 * @param head existing BlastSeqLoc to append onto, if *head
 *   is NULL then it will be set to new BlastSeqLoc, may be NULL [in|out]
 * @param from Start of the interval [in]
 * @param to End of the interval [in]
 * @return Pointer to the allocated BlastSeqLoc structure (i.e.: tail of the
 * list).
 */
NCBI_XBLAST_EXPORT
BlastSeqLoc* BlastSeqLocNew(BlastSeqLoc** head, Int4 from, Int4 to);

/** Appends the BlastSeqLoc to the list of BlastSeqLoc-s pointed to by head.
 * @param head Pointer to the head of the linked list of BlastSeqLoc-s [in]
 * @param node Pointer to the node to be added to the list. If this is NULL,
 * this function does nothing. [in]
 * @returns pointer to the second argument to this function (i.e.: tail of the
 * list)
 */
NCBI_XBLAST_EXPORT
BlastSeqLoc* BlastSeqLocAppend(BlastSeqLoc** head, BlastSeqLoc* node);

/** Deallocate a single BlastSeqLoc structure and its contents, without
 * following its next pointer
 * @param node structure to deallocate [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
BlastSeqLoc* BlastSeqLocNodeFree(BlastSeqLoc* node);

/** Deallocate all BlastSeqLoc objects in a chain.
 * @param loc object to be freed [in]
 * @return NULL pointer returned.
 */
NCBI_XBLAST_EXPORT
BlastSeqLoc* BlastSeqLocFree(BlastSeqLoc* loc);

/** Make a deep copy of the linked list of BlastSeqLoc-s pointed to by its
 * argument
 * @param head head of the linked list [in]
 * @return NULL on NULL input or memory allocation failure, else a copy of the
 * list and its contents
 */
NCBI_XBLAST_EXPORT
BlastSeqLoc* BlastSeqLocListDup(BlastSeqLoc* head);

/** Converts reverse strand coordinates to forward strand in place.
 * @param masks BlastSeqLoc to be reversed [in|out]
 * @param query_length length of query [in]
 */
NCBI_XBLAST_EXPORT
void BlastSeqLocReverse(BlastSeqLoc* masks, Int4 query_length);

/** Go through all mask locations in one sequence and combine any that overlap,
 * deallocating the unneeded locations.
 * @param mask_loc The list of masks to be merged (in place) [in|out] 
 * @param link_value Largest gap size between locations for which they
 *                   should be linked together [in] 
*/
NCBI_XBLAST_EXPORT
void
BlastSeqLocCombine(BlastSeqLoc** mask_loc, Int4 link_value);

/** Allocate memory for a BlastMaskLoc.
 * @param total number of contexts for which SSeqLocs should be allocated 
 * (result of number of queries * number of contexts for given program) [in]
 * @return Pointer to the allocated BlastMaskLoc structure.
*/
NCBI_XBLAST_EXPORT
BlastMaskLoc* BlastMaskLocNew(Int4 total);

/** 
 * @brief Perform a deep copy of the BlastMaskLoc structure passed to this
 * function
 * 
 * @param mask_loc Source masking location structure [in]
 * 
 * @return Deep copy of its argument, or NULL if the argument was NULL or if
 * not enough memory was available
 */
NCBI_XBLAST_EXPORT
BlastMaskLoc* BlastMaskLocDup(const BlastMaskLoc* mask_loc);

/** Deallocate memory for a BlastMaskLoc structure
 * as well as the BlastSeqLoc's pointed to.
 * @param mask_loc the object to be deleted [in]
 * @return NULL pointer
 */
NCBI_XBLAST_EXPORT
BlastMaskLoc* BlastMaskLocFree(BlastMaskLoc* mask_loc);

/** Given a BlastMaskLoc with an array of lists of DNA mask locations, 
 * substitutes that array by a new array of per-protein-frame mask location 
 * lists.
 * @param mask_loc Mask locations structure. This structure can have either
 * masks for all frames in nucleotide coordinates (e.g.: the results of
 * translating protein masks to nucleotide) or a single mask per query
 * (i.e.:location NUM_FRAMES*query_index). In the latter case, this mask will 
 * be used for all frames. [in|out]
 * @param query_info Query information structure, containing contexts data [in]
 * @note This function does NOT take into consideration the strands requested
 * to be searched, which is INCONSISTENT with what the C++ API does.
 */
NCBI_XBLAST_EXPORT
Int2 BlastMaskLocDNAToProtein(BlastMaskLoc* mask_loc, 
                              const BlastQueryInfo* query_info);

/** Given a BlastMaskLoc with an array of lists of mask locations per protein
 * frame, recalculates all mask offsets in terms of the DNA sequence. 
 * @param mask_loc Mask locations structure [in|out]
 * @param query_info Query information structure, containing contexts data [in]
 */
NCBI_XBLAST_EXPORT
Int2 BlastMaskLocProteinToDNA(BlastMaskLoc* mask_loc, 
                              const BlastQueryInfo* query_info);

/** This function takes the list of mask locations (i.e., regions that 
 * should not be searched or not added to lookup table) and makes up a set 
 * of SSeqRange*'s in the concatenated sequence built from a set of queries,
 * that should be searched (that is, takes the complement). 
 * If all sequences in the query set are completely filtered, then an
 * SSeqRange is created and both of its elements (left and right) are set to 
 * -1 to indicate this. 
 * If any of the mask_loc's is NULL, an SSeqRange for the full span of the 
 * respective query sequence is created.
 * @param program_number Type of BLAST program [in]
 * @param query_info The query information structure [in]
 * @param mask_loc All mask locations [in]
 * @param complement_mask Linked list of SSeqRange*s in the concatenated 
 *                        sequence to be indexed in the lookup table . [out]
 */
NCBI_XBLAST_EXPORT
Int2 
BLAST_ComplementMaskLocations(EBlastProgramType program_number, 
   const BlastQueryInfo* query_info, const BlastMaskLoc* mask_loc, 
   BlastSeqLoc* *complement_mask);

/** Runs seg filtering functions, according to the filtering options, returns
 * BlastSeqLoc*. Should combine all SeqLocs so they are non-redundant.
 * @param program_number Type of BLAST program [in]
 * @param sequence The sequence or part of the sequence to be filtered [in]
 * @param length Length of the (sub)sequence [in]
 * @param offset Offset into the full sequence [in]
 * @param filter_options specifies how filtering is to be done [in]
 * @param seqloc_retval Resulting locations for filtered region. [out]
 * @param blast_message error messages on error [out]
 * @return zero on success
*/
NCBI_XBLAST_EXPORT
Int2
BlastSetUp_Filter(EBlastProgramType program_number, 
    Uint1* sequence, 
    Int4 length, 
    Int4 offset, 
    const SBlastFilterOptions* filter_options,
    BlastSeqLoc* *seqloc_retval,
    Blast_Message * *blast_message);


/** Does preparation for filtering and then calls BlastSetUp_Filter
 * @param query_blk sequence to be filtered [in]
 * @param query_info info on sequence to be filtered [in]
 * @param program_number one of blastn,blastp,blastx,etc. [in]
 * @param filter_options specifies how filtering is to be done [in]
 * @param filter_out resulting locations for filtered region. [out]
 * @param blast_message message that needs to be sent back to user.
*/
NCBI_XBLAST_EXPORT
Int2
BlastSetUp_GetFilteringLocations(BLAST_SequenceBlk* query_blk, 
                                 const BlastQueryInfo* query_info, 
                                 EBlastProgramType program_number, 
                                 const SBlastFilterOptions* filter_options, 
                                 BlastMaskLoc** filter_out, 
                                 Blast_Message* *blast_message);

/** Masks the letters in buffer.
 * This is a low-level routine and takes a raw buffer which it assumes
 * to be in ncbistdaa (protein) or blastna (nucleotide).
 * @param buffer the sequence to be masked (will be modified, cannot be NULL or
 * undefined behavior will result).[in|out]
 * @param length length of the sequence to be masked . [in]
 * @param is_na nucleotide if TRUE [in]
 * @param mask_loc the BlastSeqLoc to use for masking [in] 
 * @param reverse minus strand if TRUE [in]
 * @param offset how far along sequence is 1st residuse in buffer [in]
*/
NCBI_XBLAST_EXPORT
void
Blast_MaskTheResidues(Uint1 * buffer, Int4 length, Boolean is_na, 
    const BlastSeqLoc* mask_loc, Boolean reverse, Int4 offset);

/** Mask protein letters that are currently unsupported. This routine
 *  is used to make the core ignore letters within protein sequences
 *  that cannot (yet) be correctly handled
 * @param seq Protein sequence to be masked (ncbistdaa format required).
 *            Letters whose numerical value exceeds a cutoff are
 *            converted into kProtMask values [in|out]
 * @param min_invalid The first ncbistdaa value that is considered invalid.
 *            All sequence letters with numerical value >= this number
 *            are masked [in]
 */
NCBI_XBLAST_EXPORT
void
Blast_MaskUnsupportedAA(BLAST_SequenceBlk* seq, Uint1 min_invalid);

/** Masks the sequence given a BlastMaskLoc
 * @param query_blk sequence to be filtered [in]
 * @param query_info info on sequence to be filtered [in]
 * @param filter_maskloc Locations to filter [in]
 * @param program_number one of blastn,blastp,blastx,etc. [in]
*/
NCBI_XBLAST_EXPORT
void
BlastSetUp_MaskQuery(BLAST_SequenceBlk* query_blk, 
                     const BlastQueryInfo* query_info, 
                     const BlastMaskLoc *filter_maskloc, 
                     EBlastProgramType program_number);

/** Produces SBlastFilterOptions from a string that has been traditionally supported
 * in blast.
 * @param program_number Type of BLAST program [in]
 * @param instructions the string describing the filtering to be done [in]
 * @param filtering_options the structure to be filled in [out]
 * @param blast_message optional field for error messages [out]
 * @return zero on success
 */
NCBI_XBLAST_EXPORT
Int2
BlastFilteringOptionsFromString(EBlastProgramType program_number, 
                                const char* instructions, 
                                SBlastFilterOptions* *filtering_options, 
                                Blast_Message* *blast_message);

/// Convert the filtering options structure to a string
/// @param filtering_options filtering options structure, assumed to be
/// correctly filled in [in]
/// @return malloc'd string containing filtering options or NULL if there is
/// not enough memory to create the return value. Caller
/// must free() return value
/// @sa TAutoCharPtr it is recommended that the caller store the return value
/// in a TAutoCharPtr if working in C++ (defined in blast_aux.hpp)
NCBI_XBLAST_EXPORT
char*
BlastFilteringOptionsToString(const SBlastFilterOptions* filtering_options);

/** Determines whether this is a nucleotide query and whether this a minus strand or not
 *
 * @param is_na the query is nucleotide
 * @param context offset in the QueryInfo array
 * @return TRUE if this is minus strand
 */
static NCBI_INLINE Boolean BlastIsReverseStrand(Boolean is_na, Int4 context)
{
     return (is_na && ((context & 1) != 0));

}

#ifdef __cplusplus
}
#endif
#endif /* !ALGO_BLAST_CORE__BLAST_FILTER__H */
