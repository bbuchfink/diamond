/*  $Id: blast_encoding.h 118195 2008-01-24 21:22:19Z camacho $
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

/** @file blast_encoding.h
 *  Declarations of static arrays used to define some NCBI encodings to be used
 *  in a toolkit independent manner by the BLAST engine.
 */

#ifndef ALGO_BLAST_CORE___BLAST_ENCODING__H
#define ALGO_BLAST_CORE___BLAST_ENCODING__H

#include "ncbi_std.h"

#define NCBI_XBLAST_EXPORT

/** @addtogroup AlgoBlast
 *
 * @{
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Different types of sequence encodings for sequence retrieval from the 
 * BLAST database 
 */
typedef enum { 
    eBlastEncodingProtein       = 0, /**< NCBIstdaa */
    eBlastEncodingNucleotide    = 1, /**< Special encoding for preliminary 
                                       stage of BLAST: permutation of NCBI4na.
                                       A.k.a.: BLASTNA encoding
                                      */
    eBlastEncodingNcbi4na       = 2, /**< NCBI4na */
    eBlastEncodingNcbi2na       = 3, /**< NCBI2na */
    eBlastEncodingError         = 255 /**< Error value for encoding */
} EBlastEncoding;

/* Nucleotide encodings */

/** Translates between ncbi4na and blastna. The first four elements
 *	of this array match ncbi2na. */
NCBI_XBLAST_EXPORT extern const int NCBI4NA_TO_BLASTNA[];

/** Translates between blastna and ncbi4na. */
NCBI_XBLAST_EXPORT extern const int BLASTNA_TO_NCBI4NA[];

/** Translates between iupacna and blastna. */
NCBI_XBLAST_EXPORT extern const int IUPACNA_TO_BLASTNA[];

/** Translates between iupacna and ncbi4na. */
NCBI_XBLAST_EXPORT extern const int IUPACNA_TO_NCBI4NA[];

/** Translates between ncbieaa and ncbistdaa. */
NCBI_XBLAST_EXPORT extern const int AMINOACID_TO_NCBISTDAA[];

/** Translates between ncbieaa and ncbistdaa. */
NCBI_XBLAST_EXPORT extern const char NCBISTDAA_TO_AMINOACID[];

/** Translates between blastna and iupacna. */
NCBI_XBLAST_EXPORT extern const char BLASTNA_TO_IUPACNA[];

/** Translates between ncbi4na and iupacna. */
NCBI_XBLAST_EXPORT extern const char NCBI4NA_TO_IUPACNA[];

#define BLAST2NA_SIZE 4     /**< Size of compressed nucleic acid alphabet */
#define BLASTNA_SIZE 16     /**< Size of nucleic acid alphabet */
#define BLASTAA_SIZE 28     /**< Size of aminoacid alphabet */


#define BLASTNA_SEQ_CODE 99 /**< Identifies the blastna alphabet, for use in 
                                blast only. */
#define BLASTAA_SEQ_CODE 11 /**< == Seq_code_ncbistdaa */
#define NCBI4NA_SEQ_CODE 4  /**< == Seq_code_ncbi4na */	

/** Sentinel byte for protein sequences */
NCBI_XBLAST_EXPORT extern const int kProtSentinel;
/** Sentinel nibble for nucleotide sequences */
NCBI_XBLAST_EXPORT extern const int kNuclSentinel;

#ifdef __cplusplus
}
#endif

/* @} */

#endif  /* ALGO_BLAST_CORE___BLAST_ENCODING__H */
