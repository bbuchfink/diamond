/* $Id: blast_query_info.h 419604 2013-11-26 22:49:06Z rafanovi $
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
 * Author: Christiam Camacho
 *
 */

/** @file blast_query_info.h
 * Definitions and functions associated with the BlastQueryInfo structure
 */

#ifndef ALGO_BLAST_CORE__BLAST_QUERY_INFO__H
#define ALGO_BLAST_CORE__BLAST_QUERY_INFO__H

#include "ncbi_std.h"
#include "blast_def.h"
#include "blast_program.h"

#ifdef __cplusplus
extern "C" {
#endif

/** The context related information */
typedef struct BlastContextInfo {
    Int4 query_offset;      /**< Offset of this query, strand or frame in the
                               concatenated super-query. */
    Int4 query_length;      /**< Length of this query, strand or frame */
    Int8 eff_searchsp;      /**< Effective search space for this context. */
    Int4 length_adjustment; /**< Length adjustment for boundary conditions */
    Int4 query_index;       /**< Index of query (same for all frames) */
    Int1 frame;             /**< Frame number (-1, -2, -3, 0, 1, 2, or 3) */
    bool is_valid;       /**< Determine if this context is valid or not.
                              This field should be set only by the setup code
                              and read by subsequent stages of the BLAST search
                              */
} BlastContextInfo;

/** Forward declaration of SPHIQueryInfo */
struct SPHIQueryInfo;

/** The query related information 
 */
typedef struct BlastQueryInfo {
    Int4 first_context;  /**< Index of the first element of the context array */
    Int4 last_context;   /**< Index of the last element of the context array */
    int num_queries;     /**< Number of query sequences */
    BlastContextInfo * contexts; /**< Information per context */
    Uint4 max_length;    /**< Length of the longest among the concatenated
                            queries */
    struct SPHIQueryInfo* pattern_info; /**< Counts of PHI BLAST pattern
                                      occurrences, used in PHI BLAST only. */
} BlastQueryInfo;


/** Allocate memory for query information structure */
NCBI_XBLAST_EXPORT
BlastQueryInfo* BlastQueryInfoNew(EBlastProgramType program, int num_queries);

/** Deallocate memory for query information structure */
NCBI_XBLAST_EXPORT
BlastQueryInfo* BlastQueryInfoFree(BlastQueryInfo* query_info);

/** Duplicates the query information structure */
NCBI_XBLAST_EXPORT
BlastQueryInfo* BlastQueryInfoDup(const BlastQueryInfo* query_info);

/** Given a context from BLAST engine core, return the query index.
 * @param context Context saved in a BlastHSP structure [in]
 * @param program Type of BLAST program [in]
 * @return Query index in a set of queries or -1 on error
 */
NCBI_XBLAST_EXPORT
Int4 Blast_GetQueryIndexFromContext(Int4 context, EBlastProgramType program);


/** Return the query index (zero based), given the query offset
 *   in the initial HSP as the program.
 * @param query_offset Offset of the query in the initial HSP [in]
 * @param program EBlastProgramType [in]
 * @param query_info information about all the queries [in]
 * @return Query Index in a set of queries
*/
NCBI_XBLAST_EXPORT
Int4 Blast_GetQueryIndexFromQueryOffset(Int4 query_offset, EBlastProgramType program, const BlastQueryInfo* query_info);


/** Retrieve a query sequence's search space
 * @param qinfo BlastQueryInfo structure [in]
 * @param program CORE program type [in]
 * @param query_index number of the query 
 * (query_index < BlastQueryInfo::num_queries) [in]
 * @return the search space of the query sequence requested or 0 if this is not
 * set */
NCBI_XBLAST_EXPORT
Int8
BlastQueryInfoGetEffSearchSpace(const BlastQueryInfo* qinfo,
                                EBlastProgramType program,
                                Int4 query_index);

/** Set a query sequence's search space
 * @param qinfo BlastQueryInfo structure [in]
 * @param program CORE program type [in]
 * @param query_index number of the query 
 * (query_index < BlastQueryInfo::num_queries) [in]
 * @param eff_searchsp the effective search space to use [in]
 */
NCBI_XBLAST_EXPORT
void
BlastQueryInfoSetEffSearchSpace(BlastQueryInfo* qinfo,
                                EBlastProgramType program,
                                Int4 query_index,
                                Int8 eff_searchsp);

/** Obtains the sequence length for a given query in the query, without taking
 * into consideration any applicable translations 
 * @param qinfo BlastQueryInfo structure [in]
 * @param program CORE program type [in]
 * @param query_index number of the query 
 * (query_index < BlastQueryInfo::num_queries) [in]
 * @return the length of the query sequence requested
 */
NCBI_XBLAST_EXPORT
Int4 BlastQueryInfoGetQueryLength(const BlastQueryInfo* qinfo,
                                  EBlastProgramType program,
                                  Int4 query_index);

/** Create auxiliary query structures with all data corresponding
 * to a single query sequence within a concatenated set. Allocates the 
 * structures if the pointers are NULL on input; otherwise only changes the
 * contents.
 * @param one_query_info_ptr Pointer to the query information structure for a 
 *                           single query. Allocated and filled here, so the
 *                           caller of this function will be responsible for
 *                           freeing it. [out]
 * @param one_query_ptr Pointer to the query sequence block structure; allocated
 *                      here, but the contents are not allocated; it is still 
 *                      safe to free by the caller after use. [out]
 * @param query_info Query information structure containing information about a 
 *                   concatenated set. [in]
 * @param query Query sequence block corresponding to a concatenated set of 
 *              queries. [in]
 * @param query_index Which query index to create the auxiliary structures 
 *                    for? [in]
 * @return -1 if memory allocation failed; 0 on success
 */
NCBI_XBLAST_EXPORT
Int2 Blast_GetOneQueryStructs(BlastQueryInfo** one_query_info_ptr, 
                              BLAST_SequenceBlk** one_query_ptr,
                              const BlastQueryInfo* query_info, 
                              BLAST_SequenceBlk* query, Int4 query_index);

/** Search BlastContextInfo structures for the specified offset */
NCBI_XBLAST_EXPORT
Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A);

/** Get the number of bytes required for the concatenated sequence
 * buffer, given a query info structure.  The context data should
 * already be assigned.
 * @param qinfo  Query info structure. [in/out]
 * @return Number of bytes for all queries and inter-query marks.
 */
NCBI_XBLAST_EXPORT
Uint4
QueryInfo_GetSeqBufLen(const BlastQueryInfo* qinfo);


/** Copy the context query offsets to an allocated array of Int4.
 * @param info Describes the concatenated query.
 * @return Allocated array.
 */
NCBI_XBLAST_EXPORT
Int4 * ContextOffsetsToOffsetArray(const BlastQueryInfo* info);


/** Copy the context query offsets from an array of Int4, allocating
 * the context array if needed.
 * @param info Destination for the values.
 * @param new_offsets Array of values to copy from.
 * @param prog        The blast program type.
 */
NCBI_XBLAST_EXPORT
void OffsetArrayToContextOffsets(BlastQueryInfo    * info,
                                 Int4              * new_offsets,
                                 EBlastProgramType   prog);

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__BLAST_QUERY_INFO__H */
