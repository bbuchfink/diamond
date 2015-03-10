/* $Id: blast_message.c 403705 2013-06-18 11:19:40Z fongah2 $
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

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: blast_message.c 403705 2013-06-18 11:19:40Z fongah2 $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include "blast_def.h"
#include "blast_message.h"

/** Declared in blast_message.h as extern const. */
const int kBlastMessageNoContext = -1;

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
    
    retval = calloc(1, sizeof(SMessageOrigin));
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
    return NULL;
}

Blast_Message* 
Blast_MessageFree(Blast_Message* blast_msg)
{
        Blast_Message* var_msg = NULL;
        Blast_Message* next = NULL;

	if (blast_msg == NULL)
		return NULL;

        var_msg = blast_msg;
        while (var_msg)
        {
	     sfree(var_msg->message);
             var_msg->origin = SMessageOriginFree(var_msg->origin);
             next = var_msg->next;
	     sfree(var_msg);
             var_msg = next;
        }
        
	return NULL;
}

Int2 
Blast_MessageWrite(Blast_Message* *blast_msg, EBlastSeverity severity, 
                   int context, const char *message)
{
        Blast_Message* new_msg = NULL;

	if (blast_msg == NULL)
	     return 1;

	 new_msg = (Blast_Message*) calloc(1, sizeof(Blast_Message));
         if (new_msg == NULL)
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

Int2 
Blast_MessagePost(Blast_Message* blast_msg)
{
	if (blast_msg == NULL)
		return 1;

	fprintf(stderr, "%s", blast_msg->message);	/* FIXME! */

	return 0;
}

void
Blast_Perror(Blast_Message* *msg, Int2 error_code, int context)
{
    Blast_PerrorEx(msg, error_code, NULL, -1, context);
    return;
}

void 
Blast_PerrorEx(Blast_Message* *msg,
                              Int2 error_code, 
                              const char* file_name, 
                              int lineno,
                              int context)
{
    Blast_Message* new_msg = (Blast_Message*) calloc(1, sizeof(Blast_Message));
    ASSERT(msg);

    switch (error_code) {

    case BLASTERR_IDEALSTATPARAMCALC:
        new_msg->message = strdup("Failed to calculate ideal Karlin-Altschul "
                                 "parameters");
        new_msg->severity = eBlastSevError;
        new_msg->context = context;
        break;
    case BLASTERR_REDOALIGNMENTCORE_NOTSUPPORTED:
        new_msg->message = strdup("Composition based statistics or "
                                 "Smith-Waterman not supported for your "
                                 "program type");
        new_msg->severity = eBlastSevError;
        new_msg->context = context;
        break;
    case BLASTERR_INTERRUPTED:
        new_msg->message = strdup("BLAST search interrupted at user's request");
        new_msg->severity = eBlastSevInfo;
        new_msg->context = context;
        break;
    case BLASTERR_NOVALIDKARLINALTSCHUL:
        new_msg->message = strdup("Warning: Could not calculate ungapped Karlin-Altschul "
                               "parameters due to an invalid query sequence or its translation. "
                               "Please verify the query sequence(s) and/or filtering options");
        new_msg->severity = eBlastSevError;
        new_msg->context = context;
        break;

    /* Fatal errors */
    case BLASTERR_MEMORY:
        /** @todo Ideally this message would be more informative (the error code
         * already conveys this information) so that this string can be
         * displayed to the end user via the CATCH_ALL macro. If this string is
         * ever changed, please update that macro accordingly (ideally this
         * error code would be caught and would lead to a CBlastSystemException
         * being thrown with the eOutOfMemory error code) */
        new_msg->message = strdup("Out of memory");
        new_msg->severity = eBlastSevFatal;
        new_msg->context = context;
        break;
    case BLASTERR_INVALIDPARAM:
        new_msg->message = strdup("Invalid argument to function");
        new_msg->severity = eBlastSevFatal;
        new_msg->context = context;
        break;
    case BLASTERR_INVALIDQUERIES:
        new_msg->message = strdup("search cannot proceed due to errors in all "
                                 "contexts/frames of query sequences");
        new_msg->severity = eBlastSevFatal;
        new_msg->context = context;
        break;

    case BLASTERR_SEQSRC:
        new_msg->message = strdup("search cannot proceed due to errors "
                                 "retrieving sequences from databases");
        new_msg->severity = eBlastSevFatal;
        new_msg->context = context;
        break;
    /* No error, just free the structure */
    case 0:
        new_msg = Blast_MessageFree(new_msg);
        break;

    /* Unknown error */
    default:
        {
            char buf[512];
            snprintf(buf, sizeof(buf) - 1, "Unknown error code %d", error_code);
            new_msg->message = strdup(buf);
            new_msg->severity = eBlastSevError;
            new_msg->context = context;
        }
        break;
    }

    if (file_name && lineno > 0) {
        new_msg->origin = SMessageOriginNew(file_name, 
                                           (unsigned int) lineno);
    }

    if (*msg)
    {
          Blast_Message* var = *msg;
          while (var->next)
              var = var->next;
          var->next = new_msg;
    }
    else
    {
           *msg = new_msg;
    }

    return;
}
