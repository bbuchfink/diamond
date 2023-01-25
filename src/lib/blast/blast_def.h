/* $Id: blast_def.h 365511 2012-06-06 16:00:13Z fongah2 $
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
 *
 */

/** @file blast_def.h
 * Definitions used throughout BLAST
 */

#ifndef ALGO_BLAST_CORE__BLAST_DEF__H
#define ALGO_BLAST_CORE__BLAST_DEF__H

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>


/****************************** Constants *********************************/

extern const int kDustLevel;  /**< Level parameter used by dust. */
extern const int kDustWindow; /**< Window parameter used by dust. */
extern const int kDustLinker; /**< Parameter used by dust to link together close low-complexity segments. */

extern const int kSegWindow;  /**< Window that SEG examines at once. */
extern const double kSegLocut;   /**< Locut parameter for SEG. */
extern const double kSegHicut;   /**< Hicut parameter for SEG. */

/** Maximum number of HPSs to be saved in an ungapped search.
 * Value defined in blast_options.c
 */
extern const int kUngappedHSPNumMax;

/******************** Preprocessor definitions ******************************/

/** Codons are always of length 3 */
#ifndef CODON_LENGTH
#define CODON_LENGTH 3
#endif

/** For translated gapped searches, this is the default value in
 * nucleotides of longest_intron (for ungapped translated searches,
 * the default value of longest_intron is zero, which causes a legacy
 * method of HSP linking that does not use longest_intron to be
 * invoked).
 *
 * The value 122 corresponds to 40 amino acids: 40 codons * 3
 * nucleotides per codon + up to 2 frame shifts.  40 amino acids is
 * the maximum gap size in the untranslated sequence, so
 * DEFAULT_LONGEST_INTRON makes these two gap sizes equal.
 */ 
#ifndef DEFAULT_LONGEST_INTRON
#define DEFAULT_LONGEST_INTRON 122
#endif

/** Compression ratio of nucleotide bases (4 bases in 1 byte) */
#ifndef COMPRESSION_RATIO
#define COMPRESSION_RATIO 4
#endif

/** Number of frames to which we translate in translating searches */
#ifndef NUM_FRAMES
#define NUM_FRAMES 6
#endif

/** Number of frames in a nucleotide sequence */
#ifndef NUM_STRANDS
#define NUM_STRANDS 2
#endif

/** Length of the genetic code string */
#ifndef GENCODE_STRLEN
#define GENCODE_STRLEN 64
#endif

/**
 * A macro expression that returns 1, 0, -1 if a is greater than,
 * equal to or less than b, respectively.  This macro evaluates its
 * arguments more than once.
 */
#ifndef BLAST_CMP
#define BLAST_CMP(a,b) ((a)>(b) ? 1 : ((a)<(b) ? -1 : 0))
#endif

/** Safe free a pointer: belongs to a higher level header. */
#define sfree(x) __sfree((void**)(void*)&(x))


/** Implemented in blast_util.c. @sa sfree */
inline void
__sfree(void **x)
{
    free(*x);
    *x = NULL;
    return;
}

/** This symbol enables the verbose option in makeblastdb and other BLAST+
 * search command line applications, as well as the option to submit searches
 * to the test server in NCBI for remote BLAST searches 
#define _BLAST_DEBUG 1
*/

#if 0
/** Define this symbol to enable debugging APIs in the BlastSeqSrc interface to
 * allow diagnostics/debugging to be performed in the composition based
 * statistics code */
#define KAPPA_PRINT_DIAGNOSTICS 1
#endif

/********************* Structure definitions ********************************/

/** Structure holding a pair of offsets. Used for storing offsets for the
 * initial seeds. In most programs the offsets are query offset and subject 
 * offset of an initial word match. For PHI BLAST, the offsets are start and 
 * end of the pattern occurrence in subject, with no query information, 
 * because all pattern occurrences in subjects are aligned to all pattern 
 * occurrences in query.
 */
typedef union BlastOffsetPair {
    struct {
        uint32_t q_off;  /**< Query offset */
        uint32_t s_off;  /**< Subject offset */
    } qs_offsets;     /**< Query/subject offset pair */
    struct {
        uint32_t s_start;/**< Start offset of pattern in subject */
        uint32_t s_end;  /**< End offset of pattern in subject */
    } phi_offsets;    /**< Pattern offsets in subject (PHI BLAST only) */
} BlastOffsetPair;

/** A structure containing two integers, used e.g. for locations for the 
 * lookup table.
 */
typedef struct SSeqRange {
   int32_t left;  /**< left endpoint of range (zero based) */
   int32_t right;  /**< right endpoint of range (zero based) */
} SSeqRange;

/** Create a new SSeqRange structure with both fields initialized
 * @param start the start of the range [in]
 * @param stop the end of the range [in]
 */
SSeqRange SSeqRangeNew(int32_t start, int32_t stop);

/** Returns the index of the range, such that this element is the
 * first range that either contains the target or if no such range exists, the
 * index of the first range, such that the target is less than this range.
 * @pre ranges array is sorted on the starting coordinates (i.e.:
 * SSeqRange::left)
 * @param ranges array of SSeqRange structures to search [in]
 * @param num_ranges number of elements in the ranges array [in]
 * @param target element to look for [in]
 * @return the index of interest in the ranges array or -1 if the function was
 * called with invalid parameters
 */
int32_t
SSeqRangeArrayLessThanOrEqual(const SSeqRange* ranges, int32_t num_ranges,
                              int32_t target);

/** Used to hold a set of positions, mostly used for filtering. 
 * oid holds the index of the query sequence.
*/
typedef struct BlastSeqLoc {
        struct BlastSeqLoc *next;  /**< next in linked list */
        SSeqRange *ssr;            /**< location data on the sequence. */
} BlastSeqLoc;

/** Structure for keeping the query masking information */
typedef struct BlastMaskLoc {
   /** Total size of the BlastSeqLoc array below. This is always the number 
     of queries times the number of contexts. Note that in the case of 
     translated query searches, these locations must be provided in protein 
     coordinates to BLAST_MainSetUp.
     @sa BLAST_GetNumberOfContexts 
     @sa BlastMaskLocDNAToProtein
    */
   int32_t total_size;

   /** Array of masked locations. 
     Every query is allocated the number of contexts associated with the 
     program. In the case of nucleotide searches, the strand(s) to search 
     dictatate which elements of the array for a given query are filled. For 
     translated searches, this should also be the same (by design) but the 
     C toolkit API does NOT implement this, it rather fills all elements 
     for all queries with masked locations in protein coordinates (if any). 
     The C++ API does follow the convention which populates each element, only
     if so dictated by the strand(s) to search for each query.
     @sa BLAST_GetNumberOfContexts
    */
   BlastSeqLoc** seqloc_array; 
} BlastMaskLoc;

/** Define the possible subject masking types */
typedef enum ESubjectMaskingType {
    eNoSubjMasking,
    eSoftSubjMasking,
    eHardSubjMasking
} ESubjectMaskingType;

/** Structure to hold a sequence. */
typedef struct BLAST_SequenceBlk {
   uint8_t* sequence; /**< Sequence used for search (could be translation). */
   uint8_t* sequence_start; /**< Start of sequence, usually one byte before
                               sequence as that byte is a NULL sentinel byte.*/
   int32_t length;         /**< Length of sequence. */
   int16_t frame; /**< Frame of the query, needed for translated searches */
   int16_t subject_strand; /**< Strand of the subject sequence for translated searches.
                          Uses the same values as ENa_strand. */
   int32_t oid; /**< The ordinal id of the current sequence */
   bool sequence_allocated; /**< TRUE if memory has been allocated for
                                  sequence */
   bool sequence_start_allocated; /**< TRUE if memory has been allocated
                                        for sequence_start */
   uint8_t* sequence_start_nomask; /**< Query sequence without masking. */
   uint8_t* sequence_nomask; /**< Start of query sequence without masking. */
   bool nomask_allocated; /**< If false the two above are just pointers to
                                   sequence and sequence_start. */
   uint8_t* oof_sequence; /**< Mixed-frame protein representation of a
                             nucleotide sequence for out-of-frame alignment */
   bool oof_sequence_allocated; /**< TRUE if memory has been allocated
                                        for oof_sequence */
   uint8_t* compressed_nuc_seq; /**< 4-to-1 compressed version of sequence */
   uint8_t* compressed_nuc_seq_start; /**< start of compressed_nuc_seq */
   BlastMaskLoc* lcase_mask; /**< Locations to be masked from operations on 
                                this sequence: lookup table for query; 
                                scanning for subject. */
   bool lcase_mask_allocated; /**< TRUE if memory has been allocated for
                                    lcase_mask */
   int32_t chunk;  /**< Used for indexing only: the chunk number within the
                     subject sequence. */
   uint8_t *gen_code_string;  /**< for nucleotide subject sequences (tblast[nx]),
                              the genetic code used to create a translated
                              protein sequence (NULL if not applicable). This
                              field is NOT owned by this data structure, it's
                              owned by the genetic code singleton. 
                              @sa gencode_singleton.h
                              */
   /* BEGIN: Data members needed for masking subjects from a BLAST database */
   SSeqRange* seq_ranges;   /**< Ranges of the sequence to search */
   uint32_t num_seq_ranges;    /**< Number of elements in seq_ranges */
   bool seq_ranges_allocated;   /**< TRUE if memory has been allocated for
                                      seq_ranges */
   ESubjectMaskingType mask_type;          /**< type of subject masking */
   /* END: Data members needed for masking subjects from a BLAST database */

   uint8_t bases_offset; /* Bases offset in first byte for SRA seq */

} BLAST_SequenceBlk;

/** Information about a single pattern occurence in the query. */
typedef struct SPHIPatternInfo {
    int32_t offset;  /**< Starting offset of this pattern occurrence. */
    int32_t length;  /**< Length of this pattern occurrence. */
} SPHIPatternInfo;

/** In PHI BLAST, structure containing information about all pattern 
 * occurrences in query.
 */
typedef struct SPHIQueryInfo {
    int32_t num_patterns;  /**< Number of pattern occurrences in query. */
    SPHIPatternInfo *occurrences; /**< Array of pattern occurrence information
                                        structures. */
    int32_t allocated_size; /**< Allocated size of the occurrences array. */
    double probability; /**< Estimated probability of the pattern */
    char* pattern;   /**< Pattern used, saved here for formatting purposes. */
} SPHIQueryInfo;


/************************* Progress monitoring/interruptible API *************/

/** Enumeration for the stages in the BLAST search */
typedef enum EBlastStage {
    /** None specified */
    eNone               = 0x0,
    /** Preliminary stage */
    ePrelimSearch       = 0x1 << 0,
    /** Traceback stage */
    eTracebackSearch    = 0x1 << 1,
    /** Both preliminary and traceback stages */
    eBoth               = (ePrelimSearch | eTracebackSearch)
} EBlastStage;

/** Progress monitoring structure. This is updated by the engine to provided to
 * the user as an argument to the user-supplied callback function 
 * (TInterruptFnPtr). This function then can assess whether the search 
 * should proceed or exit prematurely.
 * @sa TInterruptFnPtr
 */
typedef struct SBlastProgress {
    EBlastStage stage;      /**< Stage of the BLAST search currently in
                              progress */
    void* user_data;        /**< Pointer to user-provided data */
} SBlastProgress;

/** Prototype for function pointer to determine whether the BLAST search
 * should proceed or be interrupted. If this function returns true, all 
 * processing must stop and the search must discard all interim results 
 * @note In order to avoid undue overhead, this function should not perform any
 * time consuming operations and should always return (i.e.: it should never 
 * block)
 */
typedef bool (*TInterruptFnPtr) (SBlastProgress* progress_info);

/** Allocates and initializes a new SBlastProgress structure.
 * @param user_data user-provided data (not owned by the resulting structure)
 * [in]
 * Implemented in blast_util.c 
 */
SBlastProgress* SBlastProgressNew(void* user_data);

/** Deallocates a SBlastProgress structure.
 * Implemented in blast_util.c */
SBlastProgress* SBlastProgressFree(SBlastProgress* progress_info);

/** Resets the progress structure to its original state (as if newly allocated)
 * for a fresh start without touching the user_data field */
void SBlastProgressReset(SBlastProgress* progress_info);

#endif /* !ALGO_BLAST_CORE__BLAST_DEF__H */
