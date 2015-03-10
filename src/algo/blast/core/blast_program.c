#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] =
    "$Id: blast_program.c 97157 2007-01-19 14:27:24Z madden $";
#endif /* SKIP_DOXYGEN_PROCESSING */
/* ===========================================================================
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
 * Author:  Christiam Camacho / Ilya Dondoshansky
 *
 */

/** @file blast_program.c
 * Implementation auxiliary functions to determine traits of the various BLAST
 * programs supported by core BLAST
 */
    
#include "blast_program.h"

/** Convert an arbitrary integer to true/false */
#define SAFE_CAST_INT_TO_BOOLEAN(p) (((p) != 0) ? TRUE : FALSE)

/* Classify query sequence */
Boolean Blast_QueryIsProtein(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & PROTEIN_QUERY_MASK); }

Boolean Blast_QueryIsNucleotide(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & NUCLEOTIDE_QUERY_MASK); }

Boolean Blast_QueryIsPssm(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & PSSM_QUERY_MASK); }

/* Classify subject sequence */
Boolean Blast_SubjectIsProtein(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & PROTEIN_SUBJECT_MASK); }

Boolean Blast_SubjectIsNucleotide(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & NUCLEOTIDE_SUBJECT_MASK); }

Boolean Blast_SubjectIsPssm(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & PSSM_SUBJECT_MASK); }

/* Handle translated searches */
Boolean Blast_QueryIsTranslated(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & TRANSLATED_QUERY_MASK); }

Boolean Blast_SubjectIsTranslated(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & TRANSLATED_SUBJECT_MASK); }

/* Handle special programs */
Boolean Blast_ProgramIsPsiBlast(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & PSSM_QUERY_MASK); }

Boolean Blast_ProgramIsPhiBlast(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & PATTERN_QUERY_MASK); }

Boolean Blast_ProgramIsRpsBlast(EBlastProgramType p)
{ return SAFE_CAST_INT_TO_BOOLEAN(p & PSSM_SUBJECT_MASK); }

Boolean Blast_ProgramIsValid(EBlastProgramType p)
{
    switch (p) {
    case eBlastTypeBlastp:
    case eBlastTypeBlastn:
    case eBlastTypeBlastx:
    case eBlastTypeTblastn:
    case eBlastTypeTblastx:
    case eBlastTypePsiBlast:
    case eBlastTypePsiTblastn:
    case eBlastTypeRpsBlast:
    case eBlastTypeRpsTblastn:
    case eBlastTypePhiBlastp:
    case eBlastTypePhiBlastn:
        return TRUE;
        break;
    default:
        return FALSE;
        break;
    }
}
