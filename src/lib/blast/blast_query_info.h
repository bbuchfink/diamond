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
#include <stdio.h>
#include <cstdint>

#ifdef __cplusplus
extern "C" {
#endif

/** Information about paired segments (for mapping short reads)
 */
typedef enum EMagicBlastSegmentInfo {
    eNoSegments = 0,             /**< Sequence is not part of a pair */
    fFirstSegmentFlag = 1,       /**< The first sequence of a pair */
    fLastSegmentFlag = 1 << 1,   /**< The last sequence of a pair */
    fPartialFlag = 1 << 2,       /**< The other segment is not present (did not
                                   pass quality filtering */
    eFirstSegment = fFirstSegmentFlag, /**< The first sequence of a pair with
                                            both sequences read and accepted */
    eLastSegment = fLastSegmentFlag    /** The last sequence of a pair with
                                           both sequences read and accepted */
} EMagicBlastSegmentInfo;

/** The context related information */
typedef struct BlastContextInfo {
    int query_offset;      /**< Offset of this query, strand or frame in the
                               concatenated super-query. */
    int query_length;      /**< Length of this query, strand or frame */
    int eff_searchsp;      /**< Effective search space for this context. */
    int length_adjustment; /**< Length adjustment for boundary conditions */
    int query_index;       /**< Index of query (same for all frames) */
    int frame;             /**< Frame number (-1, -2, -3, 0, 1, 2, or 3) */
    bool is_valid;       /**< Determine if this context is valid or not.
                              This field should be set only by the setup code
                              and read by subsequent stages of the BLAST search
                              */
    int segment_flags;      /**< Flags describing segments for paired reads */
} BlastContextInfo;

/** Forward declaration of SPHIQueryInfo */
struct SPHIQueryInfo;

/** The query related information
 */
typedef struct BlastQueryInfo {
    int first_context;  /**< Index of the first element of the context array */
    int last_context;   /**< Index of the last element of the context array */
    int num_queries;     /**< Number of query sequences */
    BlastContextInfo * contexts; /**< Information per context */
    unsigned max_length;    /**< Length of the longest among the concatenated
                            queries */
    int min_length;    /**< Length of the shortest among the concatenated
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
int Blast_GetQueryIndexFromContext(int context, EBlastProgramType program);


/** Return the query index (zero based), given the query offset
 *   in the initial HSP as the program.
 * @param query_offset Offset of the query in the initial HSP [in]
 * @param program EBlastProgramType [in]
 * @param query_info information about all the queries [in]
 * @return Query Index in a set of queries
*/
NCBI_XBLAST_EXPORT
int Blast_GetQueryIndexFromQueryOffset(int query_offset, EBlastProgramType program, const BlastQueryInfo* query_info);


/** Retrieve a query sequence's search space
 * @param qinfo BlastQueryInfo structure [in]
 * @param program CORE program type [in]
 * @param query_index number of the query
 * (query_index < BlastQueryInfo::num_queries) [in]
 * @return the search space of the query sequence requested or 0 if this is not
 * set */
NCBI_XBLAST_EXPORT
unsigned
BlastQueryInfoGetEffSearchSpace(const BlastQueryInfo* qinfo,
                                EBlastProgramType program,
                                int query_index);

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
                                int query_index,
                                int eff_searchsp);

/** Obtains the sequence length for a given query in the query, without taking
 * into consideration any applicable translations
 * @param qinfo BlastQueryInfo structure [in]
 * @param program CORE program type [in]
 * @param query_index number of the query
 * (query_index < BlastQueryInfo::num_queries) [in]
 * @return the length of the query sequence requested
 */
NCBI_XBLAST_EXPORT
int BlastQueryInfoGetQueryLength(const BlastQueryInfo* qinfo,
                                  EBlastProgramType program,
                                  int query_index);

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
int Blast_GetOneQueryStructs(BlastQueryInfo** one_query_info_ptr,
                              BLAST_SequenceBlk** one_query_ptr,
                              const BlastQueryInfo* query_info,
                              BLAST_SequenceBlk* query, int query_index);

/** Search BlastContextInfo structures for the specified offset */
NCBI_XBLAST_EXPORT
int BSearchContextInfo(int n, const BlastQueryInfo * A);

/** Get the number of bytes required for the concatenated sequence
 * buffer, given a query info structure.  The context data should
 * already be assigned.
 * @param qinfo  Query info structure. [in/out]
 * @return Number of bytes for all queries and inter-query marks.
 */
NCBI_XBLAST_EXPORT
int
QueryInfo_GetSeqBufLen(const BlastQueryInfo* qinfo);


/** Copy the context query offsets to an allocated array of int.
 * @param info Describes the concatenated query.
 * @return Allocated array.
 */
NCBI_XBLAST_EXPORT
int * ContextOffsetsToOffsetArray(const BlastQueryInfo* info);


/** Copy the context query offsets from an array of int, allocating
 * the context array if needed.
 * @param info Destination for the values.
 * @param new_offsets Array of values to copy from.
 * @param prog        The blast program type.
 */
NCBI_XBLAST_EXPORT
void OffsetArrayToContextOffsets(BlastQueryInfo    * info,
                                 int              * new_offsets,
                                 EBlastProgramType   prog);

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__BLAST_QUERY_INFO__H */

