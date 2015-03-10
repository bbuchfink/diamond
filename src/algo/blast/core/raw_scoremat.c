/*  $Id: raw_scoremat.c 138208 2008-08-22 14:43:58Z ucko $
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
 * Author:  Aaron Ucko
 *
 * File Description:
 *   Protein alignment score matrices; shared between the two toolkits.
 *
 */

#include "raw_scoremat.h"

#include <ctype.h>
#include <string.h>

#include "sm_blosum45.c"
#include "sm_blosum50.c"
#include "sm_blosum62.c"
#include "sm_blosum80.c"
#include "sm_blosum90.c"
#include "sm_pam30.c"
#include "sm_pam70.c"
#include "sm_pam250.c"

static const char kNCBIstdaa[] = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";


int NCBISM_GetIndex(const SNCBIPackedScoreMatrix* sm, int aa)
{
    const char *p;

    /* Translate to NCBIeaa */
    if (aa >= 0  &&  aa < sizeof(kNCBIstdaa)) {
        aa = kNCBIstdaa[aa];
    } else if (islower((unsigned char) aa)) {
        aa = toupper((unsigned char) aa);
    }

    p = strchr(sm->symbols, aa);
    return p ? p - sm->symbols : -1;
}


TNCBIScore NCBISM_GetScore(const SNCBIPackedScoreMatrix* sm,
                           int aa1, int aa2)
{
    int i1, i2;
    i1 = NCBISM_GetIndex(sm, aa1);
    i2 = NCBISM_GetIndex(sm, aa2);
    if (i1 >=0  &&  i2 >= 0) {
        return sm->scores[i1 * strlen(sm->symbols) + i2];
    } else {
        return sm->defscore;
    }
}


void NCBISM_Unpack(const SNCBIPackedScoreMatrix* psm,
                   SNCBIFullScoreMatrix* fsm)
{
    const char* sym;
    int         dim, i, j, aa1, aa2;

    sym = psm->symbols;
    dim = strlen(sym);
    /* fill with default */
    memset(&fsm->s, psm->defscore, NCBI_FSM_DIM * NCBI_FSM_DIM);
    for (i = 0;  i < dim;  ++i) {
        aa1 = sym[i];
        /* get core (NCBIeaa x NCBIeaa) */
        for (j = 0;  j < dim;  ++j) {
            aa2 = sym[j];
            fsm->s[aa1][aa2] = psm->scores[i * dim + j];
        }
        /* extend horizontally */
        for (aa2 = 0;  aa2 < sizeof(kNCBIstdaa);  ++aa2) {
            fsm->s[aa1][aa2] = fsm->s[aa1][(int)kNCBIstdaa[aa2]];
        }
        for (aa2 = 'a';  aa2 <= 'z';  ++aa2) {
            fsm->s[aa1][aa2] = fsm->s[aa1][toupper((unsigned char) aa2)];
        }
    }
    /* extend vertically */
    for (aa1 = 0;  aa1 < sizeof(kNCBIstdaa);  ++aa1) {
        memcpy(fsm->s[aa1], fsm->s[(int)kNCBIstdaa[aa1]], NCBI_FSM_DIM);
    }
    for (aa1 = 'a';  aa1 <= 'z';  ++aa1) {
        memcpy(fsm->s[aa1], fsm->s[toupper((unsigned char) aa1)], NCBI_FSM_DIM);
    }
}

static
int /* bool */ s_NCBISM_StartsWith(const char* str, const char* pfx)
{
    for ( ;  *pfx;  ++str, ++pfx) {
        if (tolower((unsigned char)*str) != *pfx) {
            return 0;
        }
    }
    return 1;
}

const SNCBIPackedScoreMatrix* NCBISM_GetStandardMatrix(const char* name)
{
    switch (name[0]) {
    case 'B': case 'b':
        if ( !s_NCBISM_StartsWith(name, "blosum") ) {
            return NULL;
        }
        switch (name[6]) {
        case '4': return strcmp(name + 6, "45") ? NULL : &NCBISM_Blosum45;
        case '5': return strcmp(name + 6, "50") ? NULL : &NCBISM_Blosum50;
        case '6': return strcmp(name + 6, "62") ? NULL : &NCBISM_Blosum62;
        case '8': return strcmp(name + 6, "80") ? NULL : &NCBISM_Blosum80;
        case '9': return strcmp(name + 6, "90") ? NULL : &NCBISM_Blosum90;
        default:  return NULL;
        }

    case 'P': case 'p':
        if ( !s_NCBISM_StartsWith(name, "pam") ) {
            return NULL;
        }
        switch (name[3]) {
        case '2': return strcmp(name + 3, "250") ? NULL : &NCBISM_Pam250;
        case '3': return strcmp(name + 3, "30")  ? NULL : &NCBISM_Pam30;
        case '7': return strcmp(name + 3, "70")  ? NULL : &NCBISM_Pam70;
        }

    default:
        return NULL;
    }
}
