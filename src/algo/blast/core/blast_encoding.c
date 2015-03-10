/* $Id: blast_encoding.c 118195 2008-01-24 21:22:19Z camacho $
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

/** @file blast_encoding.c
 * Definitions of static arrays defined in blast_encoding.h.
 * @sa blast_encoding.h
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] =
    "$Id: blast_encoding.c 118195 2008-01-24 21:22:19Z camacho $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include "blast_encoding.h"

const Uint1 NCBI4NA_TO_BLASTNA[BLASTNA_SIZE] = {
    15,     /* Gap, 0 */
     0,     /* A,   1 */
     1,     /* C,   2 */
     6,     /* M,   3 */
     2,     /* G,   4 */
     4,     /* R,   5 */
     9,     /* S,   6 */
    13,     /* V,   7 */
     3,     /* T,   8 */
     8,     /* W,   9 */
     5,     /* Y,  10 */
    12,     /* H,  11 */
     7,     /* K,  12 */
    11,     /* D,  13 */
    10,     /* B,  14 */
    14      /* N,  15 */
};

const Uint1 BLASTNA_TO_NCBI4NA[BLASTNA_SIZE] = {
     1,     /* A, 0 */
     2,     /* C, 1 */
     4,     /* G, 2 */
     8,     /* T, 3 */
     5,     /* R, 4 */
    10,     /* Y, 5 */
     3,     /* M, 6 */
    12,     /* K, 7 */
     9,     /* W, 8 */
     6,     /* S, 9 */
    14,     /* B, 10 */
    13,     /* D, 11 */
    11,     /* H, 12 */
     7,     /* V, 13 */
    15,     /* N, 14 */
     0      /* Gap, 15 */
};

const char BLASTNA_TO_IUPACNA[BLASTNA_SIZE] = {
    'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 
    'W', 'S', 'B', 'D', 'H', 'V', 'N', '-'
};

const char NCBI4NA_TO_IUPACNA[BLASTNA_SIZE] = {
    '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
    'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
};

const Uint1 IUPACNA_TO_BLASTNA[128]={
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15, 0,10, 1,11,15,15, 2,12,15,15, 7,15, 6,14,15,
15,15, 4, 9, 3,15,13, 8,15, 5,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15};

const Uint1 IUPACNA_TO_NCBI4NA[128]={
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
 0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const Uint1 AMINOACID_TO_NCBISTDAA[128] = {
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,25, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,27,10,11,12,13,26,
14,15,16,17,18,24,19,20,21,22,23, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
'-','A','B','C','D','E','F','G','H','I','K','L','M',
'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
'O', 'J'};

const Uint1 kProtSentinel = NULLB;
const Uint1 kNuclSentinel = 0xF;
