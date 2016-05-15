/*  $Id: sm_pam250.c 90506 2006-09-25 19:30:59Z madden $
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

/** Entries for the PAM250 matrix at a scale of ln(2)/3.0. */

static const TNCBIScore s_Pam250PSM[25 * 25] = {
    /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,
             F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */ 
    /*A*/    2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1,
            -3,  1,  1,  1, -6, -3,  0,  0, -1,  0, -1, -8,
    /*R*/   -2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0,
            -4,  0,  0, -1,  2, -4, -2, -1, -3,  0, -1, -8,
    /*N*/    0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2,
            -3,  0,  1,  0, -4, -2, -2,  2, -3,  1, -1, -8,
    /*D*/    0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3,
            -6, -1,  0,  0, -7, -4, -2,  3, -3,  3, -1, -8,
    /*C*/   -2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5,
            -4, -3,  0, -2, -8,  0, -2, -4, -5, -5, -1, -8,
    /*Q*/    0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1,
            -5,  0, -1, -1, -5, -4, -2,  1, -2,  3, -1, -8,
    /*E*/    0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2,
            -5, -1,  0,  0, -7, -4, -2,  3, -3,  3, -1, -8,
    /*G*/    1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3,
            -5,  0,  1,  0, -7, -5, -1,  0, -4,  0, -1, -8,
    /*H*/   -1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2,
            -2,  0, -1, -1, -3,  0, -2,  1, -2,  2, -1, -8,
    /*I*/   -1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,
             1, -2, -1,  0, -5, -1,  4, -2,  3, -2, -1, -8,
    /*L*/   -2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,
             2, -3, -3, -2, -2, -1,  2, -3,  5, -3, -1, -8,
    /*K*/   -1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0,
            -5, -1,  0,  0, -3, -4, -2,  1, -3,  0, -1, -8,
    /*M*/   -1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,
             0, -2, -2, -1, -4, -2,  2, -2,  3, -2, -1, -8,
    /*F*/   -3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,
             9, -5, -3, -3,  0,  7, -1, -4,  2, -5, -1, -8,
    /*P*/    1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2,
            -5,  6,  1,  0, -6, -5, -1, -1, -2,  0, -1, -8,
    /*S*/    1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2,
            -3,  1,  2,  1, -2, -3, -1,  0, -2,  0, -1, -8,
    /*T*/    1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1,
            -3,  0,  1,  3, -5, -3,  0,  0, -1, -1, -1, -8,
    /*W*/   -6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,
             0, -6, -2, -5, 17,  0, -6, -5, -3, -6, -1, -8,
    /*Y*/   -3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,
             7, -5, -3, -3,  0, 10, -2, -3, -1, -4, -1, -8,
    /*V*/    0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2,
            -1, -1, -1,  0, -6, -2,  4, -2,  2, -2, -1, -8,
    /*B*/    0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2,
            -4, -1,  0,  0, -5, -3, -2,  3, -3,  2, -1, -8,
    /*J*/   -1, -3, -3, -3, -5, -2, -3, -4, -2,  3,  5, -3,  3,
             2, -2, -2, -1, -3, -1,  2, -3,  5, -2, -1, -8,
    /*Z*/    0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2,
            -5,  0,  0, -1, -6, -4, -2,  2, -2,  3, -1, -8,
    /*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -8,
    /***/   -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,
            -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1
};
const SNCBIPackedScoreMatrix NCBISM_Pam250 = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_Pam250PSM,
    -8
};

