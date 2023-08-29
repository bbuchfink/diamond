/* $Id$
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
 */

/** @file blast_message.c
 * These functions provide access to Blast_Message objects, used by
 * the BLAST code as a wrapper for error and warning messages.
 */

#include "blast_def.h"
#include "blast_message.h"
#include <string.h>

/** Declared in blast_message.h as extern const. */
const int kBlastMessageNoContext = -1;
const char* kBlastErrMsg_CantCalculateUngappedKAParams
    = "Could not calculate ungapped Karlin-Altschul parameters due "
      "to an invalid query sequence or its translation. Please verify the "
      "query sequence(s) and/or filtering options";

/** Allocate a new SMessageOrigin structure
 * @param filename name of the file [in]
 * @param lineno line number in the file above [in]
 * @return newly allocated structure or NULL in case of memory allocation
 * failure.
 */
SMessageOrigin* SMessageOriginNew(const char* filename, unsigned int lineno)
{
    SMessageOrigin* retval = NULL;

    if ( !filename || !(strlen(filename) > 0) ) {
        return NULL;
    }

    if ( !retval ) {
        return NULL;
    }

    retval->filename = strdup(filename);
    retval->lineno = lineno;
    return retval;
}

/** Deallocate a SMessageOrigin structure
 * @param msgo structure to deallocate [in]
 * @return NULL
 */
SMessageOrigin* SMessageOriginFree(SMessageOrigin* msgo)
{
    if (msgo) {
        sfree(msgo->filename);
        sfree(msgo);
    }
    return 0;
}

Blast_Message* 
Blast_MessageFree(Blast_Message* blast_msg)
{
        Blast_Message* var_msg = 0;
        Blast_Message* next = 0;

	if (blast_msg == 0)
		return 0;

        var_msg = blast_msg;
        while (var_msg)
        {
	     sfree(var_msg->message);
             var_msg->origin = SMessageOriginFree(var_msg->origin);
             next = var_msg->next;
	     sfree(var_msg);
             var_msg = next;
        }
        
	return 0;
}

int
Blast_MessageWrite(Blast_Message* *blast_msg, EBlastSeverity severity, 
                   int context, const char *message)
{
        Blast_Message* new_msg = 0;

	if (blast_msg == 0)
	     return 1;

	 new_msg = (Blast_Message*) calloc(1, sizeof(Blast_Message));
         if (new_msg == 0)
             return -1;

	new_msg->severity = severity;
	new_msg->context = context;
	new_msg->message = strdup(message);

        if (*blast_msg)
        {
           Blast_Message* var_msg = *blast_msg;
           while (var_msg->next)
           {
                 var_msg = var_msg->next;
           }
           var_msg->next = new_msg;
        }
        else
        {
           *blast_msg = new_msg;
        }

	return 0;
}

int
Blast_MessagePost(Blast_Message* blast_msg)
{
	if (blast_msg == 0)
		return 1;


	return 0;
}


