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

#include "blast_def.h"

/** BLASTNA element used to mask bases in BLAST */
extern const uint8_t kNuclMask;
/** NCBISTDAA element used to mask residues in BLAST */
extern const uint8_t kProtMask;

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
BlastSeqLoc* BlastSeqLocNew(BlastSeqLoc** head, uint32_t from, uint32_t to);

/** Appends the BlastSeqLoc to the list of BlastSeqLoc-s pointed to by head.
 * @param head Pointer to the head of the linked list of BlastSeqLoc-s [in]
 * @param node Pointer to the node to be added to the list. If this is NULL,
 * this function does nothing. [in]
 * @returns pointer to the second argument to this function (i.e.: tail of the
 * list)
 */
BlastSeqLoc* BlastSeqLocAppend(BlastSeqLoc** head, BlastSeqLoc* node);

/** Deallocate a single BlastSeqLoc structure and its contents, without
 * following its next pointer
 * @param node structure to deallocate [in]
 * @return NULL
 */
BlastSeqLoc* BlastSeqLocNodeFree(BlastSeqLoc* node);

/** Deallocate all BlastSeqLoc objects in a chain.
 * @param loc object to be freed [in]
 * @return NULL pointer returned.
 */
BlastSeqLoc* BlastSeqLocFree(BlastSeqLoc* loc);

/** Make a deep copy of the linked list of BlastSeqLoc-s pointed to by its
 * argument
 * @param head head of the linked list [in]
 * @return NULL on NULL input or memory allocation failure, else a copy of the
 * list and its contents
 */
BlastSeqLoc* BlastSeqLocListDup(BlastSeqLoc* head);

/** Converts reverse strand coordinates to forward strand in place.
 * @param masks BlastSeqLoc to be reversed [in|out]
 * @param query_length length of query [in]
 */
void BlastSeqLocReverse(BlastSeqLoc* masks, uint32_t query_length);

/** Go through all mask locations in one sequence and combine any that overlap,
 * deallocating the unneeded locations.
 * @param mask_loc The list of masks to be merged (in place) [in|out] 
 * @param link_value Largest gap size between locations for which they
 *                   should be linked together [in] 
*/
void
BlastSeqLocCombine(BlastSeqLoc** mask_loc, uint32_t link_value);

/** Allocate memory for a BlastMaskLoc.
 * @param total number of contexts for which SSeqLocs should be allocated 
 * (result of number of queries * number of contexts for given program) [in]
 * @return Pointer to the allocated BlastMaskLoc structure.
*/
BlastMaskLoc* BlastMaskLocNew(uint32_t total);

/** 
 * @brief Perform a deep copy of the BlastMaskLoc structure passed to this
 * function
 * 
 * @param mask_loc Source masking location structure [in]
 * 
 * @return Deep copy of its argument, or NULL if the argument was NULL or if
 * not enough memory was available
 */
BlastMaskLoc* BlastMaskLocDup(const BlastMaskLoc* mask_loc);

/** Deallocate memory for a BlastMaskLoc structure
 * as well as the BlastSeqLoc's pointed to.
 * @param mask_loc the object to be deleted [in]
 * @return NULL pointer
 */
BlastMaskLoc* BlastMaskLocFree(BlastMaskLoc* mask_loc);

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
void
Blast_MaskTheResidues(uint8_t * buffer, uint32_t length, bool is_na,
    const BlastSeqLoc* mask_loc, bool reverse, uint32_t offset);

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
void
Blast_MaskUnsupportedAA(BLAST_SequenceBlk* seq, uint8_t min_invalid);

#endif /* !ALGO_BLAST_CORE__BLAST_FILTER__H */
