/* $Id: blast_query_info.c 419604 2013-11-26 22:49:06Z rafanovi $
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

/** @file blast_query_info.c
 * Functions to manipulate the BlastQueryInfo structure
 */


#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: blast_query_info.c 419604 2013-11-26 22:49:06Z rafanovi $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include "blast_util.h"
#include "blast_query_info.h"

Int4 
Blast_GetQueryIndexFromContext(Int4 context, EBlastProgramType program)
{
   if (program == eBlastTypePsiTblastn || Blast_QueryIsProtein(program)) {
           return context;
   } else if (Blast_QueryIsTranslated(program)) {
           return context / NUM_FRAMES;
   } else {
           return context / NUM_STRANDS;
   }
}

Int4
Blast_GetQueryIndexFromQueryOffset(Int4 query_offset, EBlastProgramType program, const BlastQueryInfo* query_info)
{
   int context = BSearchContextInfo(query_offset, query_info);
  
   return Blast_GetQueryIndexFromContext(context, program);
}

BlastQueryInfo* BlastQueryInfoNew(EBlastProgramType program, int num_queries)
{
    const unsigned int kNumContexts = BLAST_GetNumberOfContexts(program);
    BlastQueryInfo* retval = NULL;
    
    if (num_queries <= 0) {
        return retval;
    }
    ASSERT(kNumContexts != 0);

    retval = (BlastQueryInfo*) calloc(1, sizeof(BlastQueryInfo));
    if ( !retval ) {
        return BlastQueryInfoFree(retval);
    }

    retval->num_queries = num_queries;

    retval->first_context = 0;
    retval->last_context = retval->num_queries * kNumContexts - 1;

    retval->contexts = (BlastContextInfo*) calloc(retval->last_context + 1,
                                                  sizeof(BlastContextInfo));

    if ( !retval->contexts ) {
        return BlastQueryInfoFree(retval);
    } else {
        int i;
        for (i = 0; i < retval->last_context + 1; i++) {
            retval->contexts[i].query_index =
                Blast_GetQueryIndexFromContext(i, program);
            ASSERT(retval->contexts[i].query_index != -1);

            retval->contexts[i].frame = BLAST_ContextToFrame(program,  i);
            ASSERT(retval->contexts[i].frame != INT1_MAX);

            retval->contexts[i].is_valid = TRUE;
        }
    }

    return retval;
}

BlastQueryInfo* BlastQueryInfoFree(BlastQueryInfo* query_info)
{
    if (query_info) {
        sfree(query_info->contexts);
        assert(query_info->pattern_info == 0);
        sfree(query_info);
    }
    return NULL;
}

BlastQueryInfo* BlastQueryInfoDup(const BlastQueryInfo* query_info)
{
   BlastQueryInfo* retval = BlastMemDup(query_info, sizeof(BlastQueryInfo));
   Int4 num_contexts = query_info->last_context + 1;
   
   retval->contexts =
       BlastMemDup(query_info->contexts, num_contexts * sizeof(BlastContextInfo));
   
   if (query_info->pattern_info) {
       assert(retval->pattern_info == 0);
   }

   return retval;
}

/** Calculates length of the DNA query from the BlastQueryInfo structure that 
 * contains context information for translated frames for a set of queries.
 * @param query_info Query information containing data for all contexts [in]
 * @param query_index Which query to find DNA length for?
 * @return DNA length of the query, calculated as sum of 3 protein frame lengths, 
 *         plus 2, because 2 last nucleotide residues do not have a 
 *         corresponding codon.
 */
static Int4 
s_GetTranslatedQueryDNALength(const BlastQueryInfo* query_info, Int4 query_index)
{
    Int4 start_context = NUM_FRAMES*query_index;
    Int4 dna_length = 2;
    Int4 index;
 
    /* Make sure that query index is within appropriate range, and that this is
       really a translated search */
    ASSERT(query_index < query_info->num_queries);
    ASSERT(start_context < query_info->last_context);

    /* If only reverse strand is searched, then forward strand contexts don't 
       have lengths information */
    if (query_info->contexts[start_context].query_length == 0)
        start_context += 3;

    for (index = start_context; index < start_context + 3; ++index)
        dna_length += query_info->contexts[index].query_length;
 
    return dna_length;
}

Int4 BlastQueryInfoGetQueryLength(const BlastQueryInfo* qinfo,
                                  EBlastProgramType program,
                                  Int4 query_index)
{
    const Uint4 kNumContexts = BLAST_GetNumberOfContexts(program);
    ASSERT(query_index < qinfo->num_queries);

    if (Blast_QueryIsTranslated(program)) {
        return s_GetTranslatedQueryDNALength(qinfo, query_index);
    } else if (program == eBlastTypeBlastn) {
        Int4 retval = qinfo->contexts[query_index*kNumContexts].query_length;
        if (retval <= 0) {
            retval = qinfo->contexts[query_index*kNumContexts+1].query_length;
        }
        return retval;
    } else {
        return qinfo->contexts[query_index*kNumContexts].query_length;
    }
}

/* FIXME: should the EBlastProgramType be added as a member of the
 * BlastQueryInfo structure? Without it, there's many operations that can't be
 * done, so it doesn't make sense to have them separate... */
Int8
BlastQueryInfoGetEffSearchSpace(const BlastQueryInfo* qinfo,
                                EBlastProgramType program,
                                Int4 query_index)
{
    Int8 retval = 0;
    Int4 i = 0;
    const Int4 kNumContexts = (Int4)BLAST_GetNumberOfContexts(program);
    ASSERT(query_index < qinfo->num_queries);

    for (i = query_index*kNumContexts; i < (query_index+1)*kNumContexts; i++) {
        if ( (retval = qinfo->contexts[i].eff_searchsp) != 0) {
            break;
        }
    }
    return retval;
}

void
BlastQueryInfoSetEffSearchSpace(BlastQueryInfo* qinfo,
                                EBlastProgramType program,
                                Int4 query_index,
                                Int8 eff_searchsp)
{
    Int4 i = 0;
    const Int4 kNumContexts = (Int4)BLAST_GetNumberOfContexts(program);
    ASSERT(query_index < qinfo->num_queries);

    for (i = query_index*kNumContexts; i < (query_index+1)*kNumContexts; i++) {
        qinfo->contexts[i].eff_searchsp = eff_searchsp;
    }
}

Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
{
    Int4 m=0, b=0, e=0, size=0;
    
    size = A->last_context+1;
    
    b = 0;
    e = size;
    while (b < e - 1) {
	m = (b + e) / 2;
	if (A->contexts[m].query_offset > n)
	    e = m;
	else
	    b = m;
    }
    return b;
}

Uint4
QueryInfo_GetSeqBufLen(const BlastQueryInfo* qinfo)
{
    BlastContextInfo * cinfo = & qinfo->contexts[qinfo->last_context];
    return cinfo->query_offset + cinfo->query_length + (cinfo->query_length ? 2 : 1);
}

Int4 *
ContextOffsetsToOffsetArray(const BlastQueryInfo* info)
{
    /* The Many Values of 'Length'
     *
     * 1. info->last_context: the index of the last query offset.
     *
     * 2. count: the number of query offsets.
     *
     * 3. count + 1: the size of the output array (has an 'extra'
     *    member so as to communicate the last sequence length).
     *
     * 4. sz: the size of the context object array
     */
    
    Uint4   count  = (info->last_context + 1);
    Uint4   sz     = sizeof(Int4) * (count+1);
    Uint4   frame  = 0;
    Int4  * result = 0;
    
    ASSERT(info);
    ASSERT(info->contexts);
    
    result = malloc(sz);
    memset(result, 0, sz);
    
    for(frame = 0; frame < count; frame++) {
        result[frame] = info->contexts[frame].query_offset;
    }
    
    /* One more entry, provides length info for last element. */
    
    result[count] = info->contexts[count-1].query_offset;
    
    if (info->contexts[count-1].query_length) {
        result[count] += info->contexts[count-1].query_length + 1;
    }
    
    return result;
}

void
OffsetArrayToContextOffsets(BlastQueryInfo    * info,
                            Int4              * new_offsets,
                            EBlastProgramType   prog)
{
    Uint4 count = (info->last_context + 1);
    Uint4 i     = 0;
    
    ASSERT(info);
    ASSERT(new_offsets);
    
    if (! info->contexts) {
        info->contexts = calloc(count, sizeof(BlastContextInfo));
    }
    
    for(i = 0; i < count; i++) {
        Int4 distance = 0;
        
        info->contexts[i].query_offset = new_offsets[i];
        
        distance = new_offsets[i+1] - new_offsets[i];
        info->contexts[i].query_length = distance ? distance-1 : 0;
        
        /* Set the frame and query index */
        
        info->contexts[i].frame =
            BLAST_ContextToFrame(prog, i);
        
        info->contexts[i].query_index =
            Blast_GetQueryIndexFromContext(i, prog);
    }
}

Int2
Blast_GetOneQueryStructs(BlastQueryInfo** one_query_info_ptr, 
                         BLAST_SequenceBlk** one_query_ptr,
                         const BlastQueryInfo* query_info, 
                         BLAST_SequenceBlk* query, Int4 query_index)
{
    Int4 num_frames;
    Int4 index;
    Int4 first_context;
    Int4 query_offset;
    BlastQueryInfo* one_query_info = NULL;
    BLAST_SequenceBlk* one_query = NULL;

    if (!one_query_info_ptr || !one_query_ptr || !query_info || !query ||
        query_index >= query_info->num_queries)
        return -1;

    num_frames = (query_info->last_context / query_info->num_queries) + 1;
    first_context = query_index*num_frames;
    query_offset = query_info->contexts[first_context].query_offset;

    one_query_info = *one_query_info_ptr;
    /* If this hasn't been already done, allocate new query information 
       structure. */
    if (!one_query_info) {
        one_query_info = (BlastQueryInfo*) calloc(1, sizeof(BlastQueryInfo));
        *one_query_info_ptr = one_query_info;
        one_query_info->contexts = (BlastContextInfo*) calloc(num_frames, sizeof(BlastContextInfo));
    }
    one_query = *one_query_ptr;
    /* If this hasn't been already done, allocate new sequence block. */
    if (!one_query) {
        one_query = (BLAST_SequenceBlk*) calloc(1, sizeof(BLAST_SequenceBlk));
        *one_query_ptr = one_query;
    }
    if (!one_query_info || !one_query)
        return -1;

    one_query_info->num_queries = 1;
    one_query_info->last_context = num_frames - 1;
    
    memcpy(one_query_info->contexts,
           &query_info->contexts[first_context],
           num_frames * sizeof(BlastContextInfo));
    
    /* Make context offsets relative to this query. */
    for (index = 0; index < num_frames; ++index) {
        one_query_info->contexts[index].query_offset -= query_offset;
    }
    
    /* Fill the sequence block information for this one query. */
    memset(one_query, 0, sizeof(BLAST_SequenceBlk));
    one_query->sequence = &query->sequence[query_offset];
    one_query->length =
        one_query_info->contexts[num_frames-1].query_offset +
        one_query_info->contexts[num_frames-1].query_length;
    one_query->sequence_allocated = FALSE;
    one_query->oid = query_index;

    return 0;
}

