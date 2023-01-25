/* $Id: ncbi_std.h 417725 2013-11-08 18:44:19Z rafanovi $
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
 *
 */

/** @file ncbi_std.h
 * Type and macro definitions from C toolkit that are not defined in C++ 
 * toolkit.
 */



#ifndef ALGO_BLAST_CORE__NCBI_STD
#define ALGO_BLAST_CORE__NCBI_STD

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <ctype.h>
#include <assert.h>

#define NCBI_XBLAST_EXPORT

/* which toolkit are we using? */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NCBI_RESTRICT /* now defined in the C++ Toolkit */
/** For some reason, ICC claims a suitable __STDC_VERSION__ but then
   barfs on restrict. Added ECC per Haruna Cofer */
#if defined(__ICC) || defined(__ECC)
#define NCBI_RESTRICT __restrict
#elif __STDC_VERSION__ >= 199901  &&  (!defined(__IBMC__) || defined(__C99_RESTRICT))
#define NCBI_RESTRICT restrict
#else
#define NCBI_RESTRICT
#endif
#endif

/* inlining support -- compiler dependent */
#if defined(__cplusplus)  ||  __STDC_VERSION__ >= 199901
/** C++ and C99 both guarantee "inline" */
#define NCBI_INLINE inline
#elif defined(__GNUC__)
/* So does GCC, normally, but it may be running with strict options
   that require the extra underscores */
#define NCBI_INLINE __inline__
#elif defined(_MSC_VER)  ||  defined(__sgi) || defined(HPUX)
/* MSVC and (older) MIPSpro always require leading underscores */
#define NCBI_INLINE __inline
#else
/** "inline" seems to work on our remaining in-house compilers
   (WorkShop, Compaq, ICC, MPW) */
#define NCBI_INLINE inline
#endif

#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strdup _strdup
#define snprintf _snprintf
#endif

#ifndef _NCBISTD_ /* if we're not in the C toolkit... */
/** bool replacment for C */
#ifndef TRUE
/** bool replacment for C indicating true. */
#define TRUE 1
#endif
#ifndef FALSE
/** bool replacment for C indicating false. */
#define FALSE 0
#endif
#endif

/** macro for assert. */
#ifndef ASSERT
#define ASSERT assert
#endif

#ifndef MIN
/** returns smaller of a and b. */
#define MIN(a,b)	((a)>(b)?(b):(a))
#endif

#ifndef MAX
/** returns larger of a and b. */
#define MAX(a,b)	((a)>=(b)?(a):(b))
#endif

#ifndef ABS
/** returns absolute value of a (|a|) */
#define ABS(a)	((a)>=0?(a):-(a))
#endif

#ifndef SIGN
/** return +1 for a > 0, -1 for a < 0 */
#define SIGN(a)	((a)>0?1:((a)<0?-1:0))
#endif

/* low-level ANSI-style functions */

#ifndef _NCBISTD_ /* if we're not in the C toolkit ... */

#ifndef UINT4_MAX
/** largest number represented by unsigned int. */
#define UINT4_MAX     4294967295U
#endif

#ifndef INT4_MAX
/** largest nubmer represented by signed int */
#define INT4_MAX    2147483647
#endif 

#ifndef INT4_MIN
/** Smallest (most negative) number represented by signed int */
#define INT4_MIN    (-2147483647-1)
#endif

#ifndef NCBIMATH_LN2
/** natural log of 2. */
#define NCBIMATH_LN2      0.69314718055994530941723212145818
#endif

#ifndef INT2_MAX
/** largest number represented by signed (two byte) short */
#define INT2_MAX    32767
#endif

#ifndef INT2_MIN
/** smallest (most negative) number represented by signed (two byte) short */
#define INT2_MIN    (-32768)
#endif

#ifndef INT1_MAX
/** largest number represented by signed short (one byte) */
#define INT1_MAX    127
#endif

#ifndef INT1_MIN
/** smallest (most negative) number represented by signed short (one byte) */
#define INT1_MIN    (-128)
#endif

#ifndef DIM
/** dimension of an array. */
#define DIM(A) (sizeof(A)/sizeof((A)[0]))
#endif

#ifndef NULLB
/** terminating byte of a char* string. */
#define NULLB '\0'
#endif

#endif /* _NCBISTD_ */

/** 64-bit integers */
#ifndef NCBI_CONST_INT8 /* C Toolkit */
#  ifdef INT64_C /* stdint.h should have this */
#    define NCBI_CONST_INT8(v)   INT64_C(v)
#    define NCBI_CONST_UINT8(v)  UINT64_C(v)
#  elif defined(_MSC_VER)
#    define NCBI_CONST_INT8(v)   v##i64
#    define NCBI_CONST_UINT8(v)  v##ui64
#  else /* Try treating as (unsigned) long long */
#    define NCBI_CONST_INT8(v)   v##LL
#    define NCBI_CONST_UINT8(v)  v##ULL
#  endif
#endif

/** Copies memory using memcpy and malloc
 * @param orig memory to be copied [in]
 * @param size amount to be copied [in]
 * @return pointer to newly allocated memory. NULL if orig NULL, size is zero,
 *   or allocation fails.
 */
NCBI_XBLAST_EXPORT
void* BlastMemDup (const void *orig, size_t size);


/******************************************************************************/

/** A generic linked list node structure */
typedef struct ListNode {
	int choice;   /**< to pick a choice */
	void *ptr;              /**< attached data */
	struct ListNode *next;  /**< next in linked list */
} ListNode;

/** Create a new list node  
 * @param vnp Pointer to the start of the list, may be NULL [in]
 * @return newly allocated node 
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeNew (ListNode* vnp);

/** Add a node to the list.
 * @param head Pointer to the start of the list, if *head is NULL will
 *  be Pointer to new node. [in] [out]
 * @return New node
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeAdd (ListNode** head);

/** Add a node to the list with a given choice and data pointer.
 * @param head Pointer to the start of the list, if *head is NULL will
 *  be Pointer to new node. [in] [out]
 * @param choice Choice value for the new node. [in]
 * @param value Data pointer for the new node. [in]
 * @return New node
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeAddPointer (ListNode** head, int choice, void *value);

/** Free all list's nodes, does not attempt to free data. 
 * @param vnp objects to be freed [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeFree (ListNode* vnp);

/** Free nodes as well as data (vnp->ptr) assuming it is one contiguous chunk.
 * @param vnp objects to be freed [in] 
 * @return NULL
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeFreeData (ListNode* vnp);

ListNode* ListNodeCopyStr (ListNode** head, int choice, const char* str);

/** Add a node to the list with a provided choice, and attached data 
 * pointing to a provided string.
 * @param head Pointer to the start of the list, if *head is NULL will
 *  be Pointer to new node. [in] [out]
 * @param choice sets the "choice" field in ListNode [in]
 * @param str char* buffer to be copied [in]
 * @return newly allocated node 
 */
NCBI_XBLAST_EXPORT

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__NCBI_STD */
