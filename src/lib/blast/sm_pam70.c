/*  $Id: sm_pam70.c 90506 2006-09-25 19:30:59Z madden $
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

/** Entries for the PAM70 matrix at a scale of ln(2)/2.0. */

static const TNCBIScore s_Pam70PSM[25 * 25] = {
    /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,
             F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */ 
    /*A*/    5, -4, -2, -1, -4, -2, -1,  0, -4, -2, -4, -4, -3,
            -6,  0,  1,  1, -9, -5, -1, -1, -3, -1, -1,-11,
    /*R*/   -4,  8, -3, -6, -5,  0, -5, -6,  0, -3, -6,  2, -2,
            -7, -2, -1, -4,  0, -7, -5, -4, -5, -2, -1,-11,
    /*N*/   -2, -3,  6,  3, -7, -1,  0, -1,  1, -3, -5,  0, -5,
            -6, -3,  1,  0, -6, -3, -5,  5, -4, -1, -1,-11,
    /*D*/   -1, -6,  3,  6, -9,  0,  3, -1, -1, -5, -8, -2, -7,
           -10, -4, -1, -2,-10, -7, -5,  5, -7,  2, -1,-11,
    /*C*/   -4, -5, -7, -9,  9, -9, -9, -6, -5, -4,-10, -9, -9,
            -8, -5, -1, -5,-11, -2, -4, -8, -7, -9, -1,-11,
    /*Q*/   -2,  0, -1,  0, -9,  7,  2, -4,  2, -5, -3, -1, -2,
            -9, -1, -3, -3, -8, -8, -4, -1, -3,  5, -1,-11,
    /*E*/   -1, -5,  0,  3, -9,  2,  6, -2, -2, -4, -6, -2, -4,
            -9, -3, -2, -3,-11, -6, -4,  2, -5,  5, -1,-11,
    /*G*/    0, -6, -1, -1, -6, -4, -2,  6, -6, -6, -7, -5, -6,
            -7, -3,  0, -3,-10, -9, -3, -1, -7, -3, -1,-11,
    /*H*/   -4,  0,  1, -1, -5,  2, -2, -6,  8, -6, -4, -3, -6,
            -4, -2, -3, -4, -5, -1, -4,  0, -4,  1, -1,-11,
    /*I*/   -2, -3, -3, -5, -4, -5, -4, -6, -6,  7,  1, -4,  1,
             0, -5, -4, -1, -9, -4,  3, -4,  4, -4, -1,-11,
    /*L*/   -4, -6, -5, -8,-10, -3, -6, -7, -4,  1,  6, -5,  2,
            -1, -5, -6, -4, -4, -4,  0, -6,  5, -4, -1,-11,
    /*K*/   -4,  2,  0, -2, -9, -1, -2, -5, -3, -4, -5,  6,  0,
            -9, -4, -2, -1, -7, -7, -6, -1, -5, -2, -1,-11,
    /*M*/   -3, -2, -5, -7, -9, -2, -4, -6, -6,  1,  2,  0, 10,
            -2, -5, -3, -2, -8, -7,  0, -6,  2, -3, -1,-11,
    /*F*/   -6, -7, -6,-10, -8, -9, -9, -7, -4,  0, -1, -9, -2,
             8, -7, -4, -6, -2,  4, -5, -7, -1, -9, -1,-11,
    /*P*/    0, -2, -3, -4, -5, -1, -3, -3, -2, -5, -5, -4, -5,
            -7,  7,  0, -2, -9, -9, -3, -4, -5, -2, -1,-11,
    /*S*/    1, -1,  1, -1, -1, -3, -2,  0, -3, -4, -6, -2, -3,
            -4,  0,  5,  2, -3, -5, -3,  0, -5, -2, -1,-11,
    /*T*/    1, -4,  0, -2, -5, -3, -3, -3, -4, -1, -4, -1, -2,
            -6, -2,  2,  6, -8, -4, -1, -1, -3, -3, -1,-11,
    /*W*/   -9,  0, -6,-10,-11, -8,-11,-10, -5, -9, -4, -7, -8,
            -2, -9, -3, -8, 13, -3,-10, -7, -5,-10, -1,-11,
    /*Y*/   -5, -7, -3, -7, -2, -8, -6, -9, -1, -4, -4, -7, -7,
             4, -9, -5, -4, -3,  9, -5, -4, -4, -7, -1,-11,
    /*V*/   -1, -5, -5, -5, -4, -4, -4, -3, -4,  3,  0, -6,  0,
            -5, -3, -3, -1,-10, -5,  6, -5,  1, -4, -1,-11,
    /*B*/   -1, -4,  5,  5, -8, -1,  2, -1,  0, -4, -6, -1, -6,
            -7, -4,  0, -1, -7, -4, -5,  5, -5,  1, -1,-11,
    /*J*/   -3, -5, -4, -7, -7, -3, -5, -7, -4,  4,  5, -5,  2,
            -1, -5, -5, -3, -5, -4,  1, -5,  5, -4, -1,-11,
    /*Z*/   -1, -2, -1,  2, -9,  5,  5, -3,  1, -4, -4, -2, -3,
            -9, -2, -2, -3,-10, -7, -4,  1, -4,  5, -1,-11,
    /*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,-11,
    /***/  -11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,
           -11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,  1
};
const SNCBIPackedScoreMatrix NCBISM_Pam70 = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_Pam70PSM,
    -11
};

