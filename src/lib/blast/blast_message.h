/* $Id: blast_message.h 403705 2013-06-18 11:19:40Z fongah2 $
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
 * Author:  Tom Madden
 *
 */

/** @file blast_message.h
 * Structures for BLAST messages
 */

#ifndef ALGO_BLAST_CORE__BLAST_MESSAGE__H
#define ALGO_BLAST_CORE__BLAST_MESSAGE__H

#include "ncbi_std.h"

#define NCBI_XBLAST_EXPORT

#ifdef __cplusplus
extern "C" {
#endif

/** Structure to enclose the origin of an error message or warning
 */
typedef struct SMessageOrigin {
    char* filename;     /**< Name of the file */
    int lineno;         /**< Line number in the file above */
} SMessageOrigin;

/** Blast error message severities .
 * These start with 1 to be consistent
 * with the C toolkit severity numbers.
 */
typedef enum {
   eBlastSevInfo = 1,
   eBlastSevWarning,
   eBlastSevError,
   eBlastSevFatal
} EBlastSeverity;

extern NCBI_XBLAST_EXPORT const int kBlastMessageNoContext;  /**< No single context is known to cause the error 
                                                 (probably a setup issue). */

/** Structure to hold the a message from the core of the BLAST engine. */
typedef struct Blast_Message {
     struct Blast_Message *next; /**< next message in this list */
     EBlastSeverity severity; /**< severity code */
     char* message;	/**< User message to be saved. */
     SMessageOrigin* origin; /**< Optional: origin of the message */
     int context; /**< Context, allows us to print message for query number. 
                      kBlastMessageNoContext used if no context applies */
} Blast_Message;

/** Deallocates message memory.
 * @param blast_msg structure to be deallocated [in]
*/

NCBI_XBLAST_EXPORT
Blast_Message* Blast_MessageFree(Blast_Message* blast_msg);


/** Writes a message to a structure.  The Blast_Message* is allocated.
 * @param blast_msg structure to be filled in [in] 
 * @param severity severity code [in] 
 * @param context query context to which this error applies [in]
 * @param message User message to be saved [in]
*/

NCBI_XBLAST_EXPORT
int Blast_MessageWrite(Blast_Message* *blast_msg, EBlastSeverity severity,
                        int context, const char *message);


/** Print a message with ErrPostEx
 * @param blast_msg message to be printed [in]
*/

NCBI_XBLAST_EXPORT
int Blast_MessagePost(Blast_Message* blast_msg);

/* FIXME: should the code below and its implementation be moved to another
 * file, say blast_error.[hc]? */

/** Analogous to perror
 * @param msg object to be appended to or created [in|out]
 * @param error_code error code returned from BLAST function [in]
 * @param context context number so that query or frame can be found [in]
 * @return Blast_Message structure containing error description
 */
NCBI_XBLAST_EXPORT
void Blast_Perror(Blast_Message* *msg, int error_code, int context);

/** Convenient define to call the function Blast_PerrorEx. */
#define Blast_PerrorWithLocation(msg, error_code, context) \
Blast_PerrorEx(msg, error_code, __FILE__, __LINE__, context)

/** Extended version of Blast_Perror which includes parameters for the file
 * name and line number where the error/warning occurred. This function should
 * be invoked via the Blast_PerrorWithLocation macro.
 * @param msg object to be appended to or created [in|out]
 * @param error_code one of the error codes defined below [in]
 * @param file_name name of the file where the error ocurred [in]
 * @param lineno line number where the error ocurred in the file above [in]
 * @param context context number so that query or frame can be found [in]
 */
NCBI_XBLAST_EXPORT
void Blast_PerrorEx(Blast_Message* *msg,
                              int error_code,
                              const char* file_name, 
                              int lineno,
                              int context);

/* BLAST error codes: these are meant to describe errors that can occur in the
 * core of BLAST only 
 */

/** System error: out of memory condition */
#define BLASTERR_MEMORY                             50

/** Invalid parameter: possible programmer error or pre-condition not met */
#define BLASTERR_INVALIDPARAM                       75

/** Could not compute the ideal Karlin-Altschul parameters */
#define BLASTERR_IDEALSTATPARAMCALC                 100

/** Composition based statistics/Smith-Waterman not supported for a program 
 * type */
#define BLASTERR_REDOALIGNMENTCORE_NOTSUPPORTED     101

/** All queries/contexts are determined invalid in the setup code */
#define BLASTERR_INVALIDQUERIES                     102

/** BLAST search was interrupted via a user-provided callback */
#define BLASTERR_INTERRUPTED                        103

/** Could not calculate Karlin-Altschul statistics for any context. */
#define BLASTERR_NOVALIDKARLINALTSCHUL              104

/** The option is not supported with the specified program. */
#define BLASTERR_OPTION_PROGRAM_INVALID             201  

  /** The value of the option is not supported (e.g., word size too small) */
#define BLASTERR_OPTION_VALUE_INVALID               202

/** Blast seqsrc returns  BLAST_SEQSRC_ERROR */
#define BLASTERR_SEQSRC								300

#ifdef __cplusplus
}
#endif
#endif /* !ALGO_BLAST_CORE__BLAST_MESSAGE__H */

