#ifndef BLAST_EXPORT__H
#define BLAST_EXPORT__H

/*  $Id$
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
 * Author:  Viatcheslav Gorelenkov
 *
 */

/** @file blast_export.h
 * Defines to provide correct exporting from BLAST DLL in Windows.
 * These are necessary to compile DLLs with Visual C++ - exports must be
 * explicitly labeled as such.
 */




#if defined(WIN32)  &&  defined(NCBI_DLL_BUILD)

#ifndef _MSC_VER
#  error "This toolkit is not buildable with a compiler other than MSVC."
#endif


#ifdef NCBI_XALGO_EXPORTS
#  define NCBI_XBLAST_EXPORT        __declspec(dllexport)
#else
#  define NCBI_XBLAST_EXPORT        __declspec(dllimport)
#endif

#elif defined(HAVE_ATTRIBUTE_VISIBILITY_DEFAULT)

#  define NCBI_XBLAST_EXPORT        __attribute__((visibility("default")))

#else

/**
 * NULL operations for other cases
 */

#  define NCBI_XBLAST_EXPORT


#endif

#endif  /*  BLAST_EXPORT__H  */
