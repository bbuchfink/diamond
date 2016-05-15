/*  $Id: sm_pam30.c 90506 2006-09-25 19:30:59Z madden $
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

/** Entries for the PAM30 matrix at a scale of ln(2)/2.0. */

static const TNCBIScore s_Pam30PSM[25 * 25] = {
    /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,
             F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */ 
    /*A*/    6, -7, -4, -3, -6, -4, -2, -2, -7, -5, -6, -7, -5,
            -8, -2,  0, -1,-13, -8, -2, -3, -6, -3, -1,-17,
    /*R*/   -7,  8, -6,-10, -8, -2, -9, -9, -2, -5, -8,  0, -4,
            -9, -4, -3, -6, -2,-10, -8, -7, -7, -4, -1,-17,
    /*N*/   -4, -6,  8,  2,-11, -3, -2, -3,  0, -5, -7, -1, -9,
            -9, -6,  0, -2, -8, -4, -8,  6, -6, -3, -1,-17,
    /*D*/   -3,-10,  2,  8,-14, -2,  2, -3, -4, -7,-12, -4,-11,
           -15, -8, -4, -5,-15,-11, -8,  6,-10,  1, -1,-17,
    /*C*/   -6, -8,-11,-14, 10,-14,-14, -9, -7, -6,-15,-14,-13,
           -13, -8, -3, -8,-15, -4, -6,-12, -9,-14, -1,-17,
    /*Q*/   -4, -2, -3, -2,-14,  8,  1, -7,  1, -8, -5, -3, -4,
           -13, -3, -5, -5,-13,-12, -7, -3, -5,  6, -1,-17,
    /*E*/   -2, -9, -2,  2,-14,  1,  8, -4, -5, -5, -9, -4, -7,
           -14, -5, -4, -6,-17, -8, -6,  1, -7,  6, -1,-17,
    /*G*/   -2, -9, -3, -3, -9, -7, -4,  6, -9,-11,-10, -7, -8,
            -9, -6, -2, -6,-15,-14, -5, -3,-10, -5, -1,-17,
    /*H*/   -7, -2,  0, -4, -7,  1, -5, -9,  9, -9, -6, -6,-10,
            -6, -4, -6, -7, -7, -3, -6, -1, -7, -1, -1,-17,
    /*I*/   -5, -5, -5, -7, -6, -8, -5,-11, -9,  8, -1, -6, -1,
            -2, -8, -7, -2,-14, -6,  2, -6,  5, -6, -1,-17,
    /*L*/   -6, -8, -7,-12,-15, -5, -9,-10, -6, -1,  7, -8,  1,
            -3, -7, -8, -7, -6, -7, -2, -9,  6, -7, -1,-17,
    /*K*/   -7,  0, -1, -4,-14, -3, -4, -7, -6, -6, -8,  7, -2,
           -14, -6, -4, -3,-12, -9, -9, -2, -7, -4, -1,-17,
    /*M*/   -5, -4, -9,-11,-13, -4, -7, -8,-10, -1,  1, -2, 11,
            -4, -8, -5, -4,-13,-11, -1,-10,  0, -5, -1,-17,
    /*F*/   -8, -9, -9,-15,-13,-13,-14, -9, -6, -2, -3,-14, -4,
             9,-10, -6, -9, -4,  2, -8,-10, -2,-13, -1,-17,
    /*P*/   -2, -4, -6, -8, -8, -3, -5, -6, -4, -8, -7, -6, -8,
           -10,  8, -2, -4,-14,-13, -6, -7, -7, -4, -1,-17,
    /*S*/    0, -3,  0, -4, -3, -5, -4, -2, -6, -7, -8, -4, -5,
            -6, -2,  6,  0, -5, -7, -6, -1, -8, -5, -1,-17,
    /*T*/   -1, -6, -2, -5, -8, -5, -6, -6, -7, -2, -7, -3, -4,
            -9, -4,  0,  7,-13, -6, -3, -3, -5, -6, -1,-17,
    /*W*/  -13, -2, -8,-15,-15,-13,-17,-15, -7,-14, -6,-12,-13,
            -4,-14, -5,-13, 13, -5,-15,-10, -7,-14, -1,-17,
    /*Y*/   -8,-10, -4,-11, -4,-12, -8,-14, -3, -6, -7, -9,-11,
             2,-13, -7, -6, -5, 10, -7, -6, -7, -9, -1,-17,
    /*V*/   -2, -8, -8, -8, -6, -7, -6, -5, -6,  2, -2, -9, -1,
            -8, -6, -6, -3,-15, -7,  7, -8,  0, -6, -1,-17,
    /*B*/   -3, -7,  6,  6,-12, -3,  1, -3, -1, -6, -9, -2,-10,
           -10, -7, -1, -3,-10, -6, -8,  6, -8,  0, -1,-17,
    /*J*/   -6, -7, -6,-10, -9, -5, -7,-10, -7,  5,  6, -7,  0,
            -2, -7, -8, -5, -7, -7,  0, -8,  6, -6, -1,-17,
    /*Z*/   -3, -4, -3,  1,-14,  6,  6, -5, -1, -6, -7, -4, -5,
           -13, -4, -5, -6,-14, -9, -6,  0, -6,  6, -1,-17,
    /*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,-17,
    /***/  -17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,
           -17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,  1
};
const SNCBIPackedScoreMatrix NCBISM_Pam30 = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_Pam30PSM,
    -17
};

