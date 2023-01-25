/* $Id: blast_program.h 103498 2007-05-04 18:56:48Z camacho $
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
 * Author:  Christiam Camacho / Ilya Dondoshansky
 *
 */

/** @file blast_program.h
 * Definitions for various programs supported by core BLAST
 */

#ifndef ALGO_BLAST_CORE___BLAST_PROGRAM__H
#define ALGO_BLAST_CORE___BLAST_PROGRAM__H

#include "ncbi_std.h"
#include <stdbool.h>


#define NCBI_XBLAST_EXPORT

#ifdef __cplusplus
extern "C" {
#endif

/* Low-level preprocessor definitions: these are not intended to be used
 * directly, use the EBlastProgramType enumeration below */

/** This bit is on if the query is protein */
#define PROTEIN_QUERY_MASK      (0x1<<0)
/** This bit is on if the subject is protein */
#define PROTEIN_SUBJECT_MASK    (0x1<<1)
/** This bit is on if the query is nucleotide */
#define NUCLEOTIDE_QUERY_MASK   (0x1<<2)
/** This bit is on if the subject is nucleotide */
#define NUCLEOTIDE_SUBJECT_MASK (0x1<<3)
/** This bit is on if the query is translated */
#define TRANSLATED_QUERY_MASK   (0x1<<4)
/** This bit is on if the subject is translated */
#define TRANSLATED_SUBJECT_MASK (0x1<<5)
/** This bit is on if the query is a PSSM (PSI-BLAST) */
#define PSSM_QUERY_MASK         (0x1<<6)
/** This bit is on if the subject is a PSSM (RPS-BLAST) */
#define PSSM_SUBJECT_MASK       (0x1<<7)
/** This bit is on if the query includes a pattern (PHI-BLAST) */
#define PATTERN_QUERY_MASK      (0x1<<8)

/******************** Main BLAST program definitions ***********************/

/** Defines the engine's notion of the different applications of the BLAST 
 * algorithm */
typedef enum {
    eBlastTypeBlastp        = (PROTEIN_QUERY_MASK | PROTEIN_SUBJECT_MASK),
    eBlastTypeBlastn        = (NUCLEOTIDE_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK),
    eBlastTypeBlastx        = (NUCLEOTIDE_QUERY_MASK | PROTEIN_SUBJECT_MASK | 
                               TRANSLATED_QUERY_MASK),
    eBlastTypeTblastn       = (PROTEIN_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK | 
                               TRANSLATED_SUBJECT_MASK),
    eBlastTypeTblastx       = (NUCLEOTIDE_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK
                               | TRANSLATED_QUERY_MASK 
                               | TRANSLATED_SUBJECT_MASK),
    eBlastTypePsiBlast      = (PSSM_QUERY_MASK | eBlastTypeBlastp),
    eBlastTypePsiTblastn    = (PSSM_QUERY_MASK | eBlastTypeTblastn),
    eBlastTypeRpsBlast      = (PSSM_SUBJECT_MASK | eBlastTypeBlastp),
    eBlastTypeRpsTblastn    = (PSSM_SUBJECT_MASK | eBlastTypeBlastx),
    eBlastTypePhiBlastp     = (PATTERN_QUERY_MASK | eBlastTypeBlastp),
    eBlastTypePhiBlastn     = (PATTERN_QUERY_MASK | eBlastTypeBlastn),
    eBlastTypeUndefined     = 0x0
} EBlastProgramType;

/************* Functions to classify programs by query **********************/

/** Returns true if the query is protein 
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_QueryIsProtein(EBlastProgramType p);

/** Returns true if the query is nucleotide
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_QueryIsNucleotide(EBlastProgramType p);

/** Returns true if the query is PSSM
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_QueryIsPssm(EBlastProgramType p);

/************* Functions to classify programs by subject **********************/

/** Returns true if the subject is protein 
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_SubjectIsProtein(EBlastProgramType p);

/** Returns true if the subject is nucleotide
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_SubjectIsNucleotide(EBlastProgramType p);

/** Returns true if the subject is PSSM
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_SubjectIsPssm(EBlastProgramType p);

/************* Functions to classify programs by translation *****************/

/** Returns true if the query is translated
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_QueryIsTranslated(EBlastProgramType p);

/** Returns true if the subject is translated
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_SubjectIsTranslated(EBlastProgramType p);

/************* Functions to classify special BLAST programs *****************/

/** Returns true if program is PSI-BLAST (i.e.: involves a PSSM as query)
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_ProgramIsPsiBlast(EBlastProgramType p);

/** Returns true if program is PHI-BLAST (i.e.: involves a pattern)
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_ProgramIsPhiBlast(EBlastProgramType p);

/** Returns true if program is RPS-BLAST (i.e.: involves a PSSM as subject)
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_ProgramIsRpsBlast(EBlastProgramType p);

/** Returns true if program is not undefined
 * @param p program type [in]
 */
NCBI_XBLAST_EXPORT 
bool Blast_ProgramIsValid(EBlastProgramType p);

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__BLAST_PROGRAM__H */

