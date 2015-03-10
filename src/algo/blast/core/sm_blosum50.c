/*  $Id: sm_blosum50.c 90507 2006-09-25 19:31:51Z madden $
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
* Author:  Aaron Ucko, Mike Gertz
*
* File Description:
*   Protein alignment score matrices; shared between the two toolkits.
*
* ===========================================================================
*/

#include "raw_scoremat.h"

/** Entries for the BLOSUM50 matrix at a scale of ln(2)/3.0. */

static const TNCBIScore s_Blosum50PSM[25 * 25] = {
    /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,
             F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */ 
    /*A*/    5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1,
            -3, -1,  1,  0, -3, -2,  0, -2, -2, -1, -1, -5,
    /*R*/   -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2,
            -3, -3, -1, -1, -3, -1, -3, -1, -3,  0, -1, -5,
    /*N*/   -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2,
            -4, -2,  1,  0, -4, -2, -3,  5, -4,  0, -1, -5,
    /*D*/   -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4,
            -5, -1,  0, -1, -5, -3, -4,  6, -4,  1, -1, -5,
    /*C*/   -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2,
            -2, -4, -1, -1, -5, -3, -1, -3, -2, -3, -1, -5,
    /*Q*/   -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0,
            -4, -1,  0, -1, -1, -1, -3,  0, -3,  4, -1, -5,
    /*E*/   -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2,
            -3, -1, -1, -1, -3, -2, -3,  1, -3,  5, -1, -5,
    /*G*/    0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3,
            -4, -2,  0, -2, -3, -3, -4, -1, -4, -2, -1, -5,
    /*H*/   -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1,
            -1, -2, -1, -2, -3,  2, -4,  0, -3,  0, -1, -5,
    /*I*/   -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,
             0, -3, -3, -1, -3, -1,  4, -4,  4, -3, -1, -5,
    /*L*/   -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,
             1, -4, -3, -1, -2, -1,  1, -4,  4, -3, -1, -5,
    /*K*/   -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2,
            -4, -1,  0, -1, -3, -2, -3,  0, -3,  1, -1, -5,
    /*M*/   -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,
             0, -3, -2, -1, -1,  0,  1, -3,  2, -1, -1, -5,
    /*F*/   -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,
             8, -4, -3, -2,  1,  4, -1, -4,  1, -4, -1, -5,
    /*P*/   -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3,
            -4, 10, -1, -1, -4, -3, -3, -2, -3, -1, -1, -5,
    /*S*/    1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2,
            -3, -1,  5,  2, -4, -2, -2,  0, -3,  0, -1, -5,
    /*T*/    0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1,
            -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -1, -5,
    /*W*/   -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,
             1, -4, -4, -3, 15,  2, -3, -5, -2, -2, -1, -5,
    /*Y*/   -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,
             4, -3, -2, -2,  2,  8, -1, -3, -1, -2, -1, -5,
    /*V*/    0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1,
            -1, -3, -2,  0, -3, -1,  5, -3,  2, -3, -1, -5,
    /*B*/   -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3,
            -4, -2,  0,  0, -5, -3, -3,  6, -4,  1, -1, -5,
    /*J*/   -2, -3, -4, -4, -2, -3, -3, -4, -3,  4,  4, -3,  2,
             1, -3, -3, -1, -2, -1,  2, -4,  4, -3, -1, -5,
    /*Z*/   -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1,
            -4, -1,  0, -1, -2, -2, -3,  1, -3,  5, -1, -5,
    /*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5,
    /***/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1
};
const SNCBIPackedScoreMatrix NCBISM_Blosum50 = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_Blosum50PSM,
    -5
};

