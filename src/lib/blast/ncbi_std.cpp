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
 */

/** @file ncbi_std.c
 * A few utilities needed by code in algo/blast/core but not provided elsewhere.
 */

#include "blast_def.h" /* for sfree() macro */
#include "ncbi_std.h"

void * BlastMemDup (const void *orig, size_t size)
{
    void*	copy;

    if (orig == NULL || size == 0)
        return NULL;

    if ((copy = malloc (size)) == NULL)
        return NULL;

    memcpy(copy, orig, size);
    return copy;
}

/*****************************************************************************
*
*   ListNodeNew(vnp)
*      adds after last node in list if vnp not NULL
*
*****************************************************************************/
ListNode* ListNodeNew (ListNode* vnp)
{
    ListNode* newnode;

    newnode = (ListNode*) calloc(1, sizeof(ListNode));
    if (vnp != NULL)
    {
        while (vnp->next != NULL)
            vnp = vnp->next;
        vnp->next = newnode;
    }
    return newnode;
}

/*****************************************************************************
*
*   ListNodeAdd(head)
*      adds after last node in list if *head not NULL
*      If *head is NULL, sets it to the new ListNode
*      returns pointer to the NEW node added
*
*****************************************************************************/
ListNode* ListNodeAdd (ListNode** head)
{
    ListNode* newnode;

    if (head != NULL)
    {
        newnode = ListNodeNew(*head);
        if (*head == NULL)
            *head = newnode;
    }
    else
        newnode = ListNodeNew(NULL);

    return newnode;
}

/*   ListNodeAddPointer (head, choice, value)
*      adds like ListNodeAdd()
*      sets newnode->choice = choice (if choice does not matter, use 0)
*      sets newnode->ptr = value
*
*****************************************************************************/
ListNode* ListNodeAddPointer (ListNode** head, unsigned choice,
                              void *value)
{
    ListNode* newnode;

    newnode = ListNodeAdd(head);
    if (newnode != NULL)
    {
        newnode->choice = choice;
        newnode->ptr = value;
    }

    return newnode;
}

/*****************************************************************************
*
*   ListNodeCopyStr (head, choice, str)
*      adds like ListNodeAdd()
*      sets newnode->choice = choice (if choice does not matter, use 0)
*      sets newnode->ptr = str
*         makes a COPY of str
*      if str == NULL, does not add a ListNode
*
*****************************************************************************/
ListNode* ListNodeCopyStr (ListNode** head, int choice, const char* str)
{
    ListNode* newnode;

    if (str == NULL) return NULL;

    newnode = ListNodeAdd(head);
    if (newnode != NULL)
    {
        newnode->choice = choice;
        newnode->ptr = strdup(str);
    }

    return newnode;
}

/*****************************************************************************
*
*   ListNodeFree(vnp)
*   	frees whole chain of ListNodes
*       Does NOT free associated data pointers
*           see ListNodeFreeData()
*
*****************************************************************************/
ListNode* ListNodeFree (ListNode* vnp)
{
    ListNode* next;

    while (vnp != NULL)
    {
        next = vnp->next;
        sfree(vnp);
        vnp = next;
    }
    return NULL;
}

/*****************************************************************************
*
*   ListNodeFreeData(vnp)
*   	frees whole chain of ListNodes
*       frees associated data pointers - BEWARE of this if these are not
*           allocated single memory block structures.
*
*****************************************************************************/
ListNode* ListNodeFreeData (ListNode* vnp)
{
    ListNode* next;

    while (vnp != NULL)
    {
        sfree(vnp->ptr);
        next = vnp->next;
        sfree(vnp);
        vnp = next;
    }
    return NULL;
}

