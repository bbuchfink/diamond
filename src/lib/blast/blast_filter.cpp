/* $Id: blast_filter.c 306966 2011-06-20 13:16:49Z maning $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
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
 */

/** @file blast_filter.c
 * All code related to query sequence masking/filtering for BLAST
 */

#include <cstddef>
#include <algorithm>
#include <assert.h>
#include "blast_filter.h"
#include "blast_seg.h"

#ifndef NULLB
/** terminating byte of a char* string. */
#define NULLB '\0'
#endif

/****************************************************************************/
/* Constants */
const uint8_t kNuclMask = 14;     /* N in BLASTNA */
const uint8_t kProtMask = 21;     /* X in NCBISTDAA */


/** Allowed length of the filtering options string. */
#define BLASTOPTIONS_BUFFER_SIZE 128

BlastSeqLoc* BlastSeqLocNew(BlastSeqLoc** head, int32_t from, int32_t to)
{
   BlastSeqLoc* loc = (BlastSeqLoc*) calloc(1, sizeof(BlastSeqLoc));
   if ( !loc ) {
       return NULL;
   }
   loc->ssr = (SSeqRange*) calloc(1, sizeof(SSeqRange));
   loc->ssr->left = from;
   loc->ssr->right = to;

   return BlastSeqLocAppend(head, loc);
}

BlastSeqLoc* BlastSeqLocAppend(BlastSeqLoc** head, BlastSeqLoc* node)
{
    if ( !node ) {
        return NULL;
    }

    if (head)
    {
        if (*head)
        {
            BlastSeqLoc* tmp = *head;
            while (tmp->next)
               tmp = tmp->next;
            tmp->next = node;
        }
        else
        {
            *head = node;
        }
    }
        
    return node;
}

/** Makes a copy of the BlastSeqLoc and also a copy of the 
 * SSRange element.  Does not copy BlastSeqLoc that is pointed
 * to by "next".
 * @param source the object to be copied [in]
 * @return another BlastSeqLoc*
 */
static BlastSeqLoc* s_BlastSeqLocNodeDup(BlastSeqLoc* source)
{
    if ( !source ) {
        return NULL;
    }
    assert(source->ssr);
    return BlastSeqLocNew(NULL, source->ssr->left, source->ssr->right);
}

/** Calculates number of links in a chain of BlastSeqLoc's.
 * @param var Chain of BlastSeqLoc structures [in]
 * @return Number of links in the chain.
 */
static int32_t s_BlastSeqLocLen(const BlastSeqLoc* var)
{
    BlastSeqLoc* itr = NULL;
    int32_t retval = 0;
   
    for (itr = (BlastSeqLoc*)var; itr; itr = itr->next, retval++) {
        ;
    }
    return retval;
}

/** Converts a BlastSeqLoc list to an array of pointers, each pointing to an
 * element of the list passed in to this function and the last element points
 * to NULL 
 * @param list List to convert to an array of pointers [in]
 * @param count number of elements populated in the array [out]
 */
static BlastSeqLoc**
s_BlastSeqLocListToArrayOfPointers(const BlastSeqLoc* list, int32_t* count)
{
    BlastSeqLoc* tmp,** retval;
    int32_t i;
    *count = 0;

    if (list == NULL) 
       return NULL;

    *count = s_BlastSeqLocLen(list);
    retval = (BlastSeqLoc**) calloc(((size_t)(*count)+1), sizeof(BlastSeqLoc*));
    for (tmp = (BlastSeqLoc*)list, i = 0; tmp != NULL && i < *count; i++) {
        retval[i] = tmp;
        tmp = tmp->next;
    }
    return retval;
}

/** Reverse elements in the list 
 * @param head pointer to pointer to the head of the list. [in|out]
 * (this is not declared static so that it can be tested in the unit tests
 */
void BlastSeqLocListReverse(BlastSeqLoc** head)
{
    BlastSeqLoc** ptrs = NULL;  /* array of pointers to BlastSeqLoc elements */
    int32_t num_elems = 0, i = 0;

    if ( !head ) {
        return;
    }

    ptrs = s_BlastSeqLocListToArrayOfPointers(*head, &num_elems);
    if (num_elems == 0) {
        return;
    }
    assert(ptrs);
    *head = ptrs[num_elems-1];
    for (i = num_elems-1; i > 0; i--) {
        ptrs[i]->next = ptrs[i-1];
    }
    ptrs[0]->next = NULL;
    sfree(ptrs);
}

BlastSeqLoc* BlastSeqLocNodeFree(BlastSeqLoc* loc)
{
    if ( !loc ) {
        return NULL;
    }
    sfree(loc->ssr);
    sfree(loc);
    return NULL;
}

BlastSeqLoc* BlastSeqLocFree(BlastSeqLoc* loc)
{
    while (loc) {
        BlastSeqLoc* next_loc = loc->next;
        loc = BlastSeqLocNodeFree(loc);
        loc = next_loc;
    }
    return NULL;
}

BlastSeqLoc* BlastSeqLocListDup(BlastSeqLoc* head)
{
    BlastSeqLoc* retval = NULL;
    BlastSeqLoc* retval_tail = NULL;

    for (; head; head = head->next) {
        retval_tail = BlastSeqLocAppend(retval_tail ? &retval_tail : &retval, 
                                        s_BlastSeqLocNodeDup(head));
    }

    return retval;
}

BlastMaskLoc* BlastMaskLocNew(int32_t total)
{
    BlastMaskLoc* retval = (BlastMaskLoc *) calloc(1, sizeof(BlastMaskLoc));
    retval->total_size = total;
    if (total > 0)
        retval->seqloc_array = (BlastSeqLoc **) calloc(total, 
                                                       sizeof(BlastSeqLoc *));
    return retval;
}

BlastMaskLoc* BlastMaskLocDup(const BlastMaskLoc* mask_loc)
{
    BlastMaskLoc* retval = NULL;
    int32_t index = 0;

    if ( !mask_loc ) {
        return NULL;
    }

    retval = BlastMaskLocNew(mask_loc->total_size);

    for (index = 0; index < mask_loc->total_size; index++) {
        retval->seqloc_array[index] =
            BlastSeqLocListDup(mask_loc->seqloc_array[index]);
    }

    return retval;
}

BlastMaskLoc* BlastMaskLocFree(BlastMaskLoc* mask_loc)
{
   int32_t index;

   if (mask_loc == NULL)
      return NULL;

   for (index=0; index<mask_loc->total_size; index++)
   {
      if (mask_loc->seqloc_array != NULL)
         BlastSeqLocFree(mask_loc->seqloc_array[index]);
   }
   sfree(mask_loc->seqloc_array);
   sfree(mask_loc);
   return NULL;
}

/** Used for qsort, compares two SeqLoc's by starting position. */
static int s_SeqRangeSortByStartPosition(const void *vp1, const void *vp2)
{
   BlastSeqLoc* v1 = *((BlastSeqLoc**) vp1);
   BlastSeqLoc* v2 = *((BlastSeqLoc**) vp2);
   SSeqRange* loc1 = (SSeqRange*) v1->ssr;
   SSeqRange* loc2 = (SSeqRange*) v2->ssr;
   
   if (loc1->left < loc2->left)
      return -1;
   else if (loc1->left > loc2->left)
      return 1;
   else
      return 0;
}

void
BlastSeqLocCombine(BlastSeqLoc** mask_loc, int32_t link_value)
{
    BlastSeqLoc** ptrs = NULL;
    int32_t i = 0, num_elems = 0;

    /* Break up the list into an array of pointers and sort it */
    ptrs = s_BlastSeqLocListToArrayOfPointers(*mask_loc, &num_elems);
    if (num_elems == 0) {
        return;
    }
    assert(ptrs);
    qsort(ptrs, (size_t)num_elems, sizeof(*ptrs), 
          s_SeqRangeSortByStartPosition);

    /* Merge the overlapping elements */
    {
        BlastSeqLoc* curr_tail = *mask_loc = ptrs[0];
        for (i = 0; i < num_elems - 1; i++) {
            const SSeqRange* next_ssr = ptrs[i+1]->ssr;
            const int32_t stop = curr_tail->ssr->right;

            if ((stop + link_value) > next_ssr->left) {
                curr_tail->ssr->right = std::max(stop, next_ssr->right);
                ptrs[i+1] = BlastSeqLocNodeFree(ptrs[i+1]);
            } else {
                curr_tail = ptrs[i+1];
            }
        }
    }

    /* Rebuild the linked list */
    {
        BlastSeqLoc* tail = *mask_loc;
        for (i = 1; i < num_elems; i++) {
            if (ptrs[i]) {
                tail->next = ptrs[i];
                tail = ptrs[i];
            }
        }
        tail->next = NULL;
    }
    sfree(ptrs);
}

void
BlastSeqLocReverse(BlastSeqLoc* masks, int32_t query_length)
{
    for(; masks; masks = masks->next) {
        masks->ssr->left    = query_length - 1 - masks->ssr->right;
        masks->ssr->right   = query_length - 1 - masks->ssr->left;
    }
}

void
Blast_MaskTheResidues(uint8_t * buffer, int32_t length, bool is_na,
                      const BlastSeqLoc* mask_loc, bool reverse, int32_t offset)
{
    const uint8_t kMaskingLetter = is_na ? kNuclMask : kProtMask;
    assert(buffer);
    for (; mask_loc; mask_loc = mask_loc->next) {

        int32_t index, start, stop;
        
        if (reverse) {
            start = length - 1 - mask_loc->ssr->right;
            stop = length - 1 - mask_loc->ssr->left;
        } else {
            start = mask_loc->ssr->left;
            stop = mask_loc->ssr->right;
        }
        
        start -= offset;
        stop -= offset;
        
        assert(start < length);
        assert(stop <= length);
        
        for (index = start; index <= stop; index++)
            buffer[index] = kMaskingLetter;
    }
}

void
Blast_MaskUnsupportedAA(BLAST_SequenceBlk* seq, uint8_t min_invalid)
{
    uint8_t *sequence = seq->sequence;
    int32_t length = seq->length;
    int32_t i;

    for (i = 0; i < length; i++) {
        if (sequence[i] >= min_invalid) {
            sequence[i] = kProtMask;
        }
    }
}
