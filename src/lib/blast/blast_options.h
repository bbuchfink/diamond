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
 *
 * Author:  Tom Madden
 *
 */

/** @file blast_options.h
 *  The structures and functions in blast_options.[ch] should be used to specify 
 *  user preferences.  The options structures should not be changed by the BLAST code
 *  but rather be read to determine user preferences.  When possible these structures
 *  should be passed in as "const".
 */

#ifndef __BLASTOPTIONS__
#define __BLASTOPTIONS__

#include "ncbi_std.h"
#include "blast_export.h"
#include "blast_program.h"
#include "blast_def.h"
#include "blast_message.h"
#include <stdbool.h>


#ifdef __cplusplus
extern "C" {
#endif


/** Some default values (used when creating blast options block and for
 * command-line program defaults. When changing these defaults, please
 * remember to update the defaults in the command-line programs 
 */

/** "window" between hits to trigger an extension. */
#define BLAST_WINDOW_SIZE_PROT 40  /**< default window (all protein searches) */
#define BLAST_WINDOW_SIZE_NUCL 0   /**< default window size (blastn) */
#define BLAST_WINDOW_SIZE_MEGABLAST 0   /**< default window size 
                                          (contiguous megablast) */
#define BLAST_WINDOW_SIZE_DISC 40  /**< default window size 
                                          (discontiguous megablast) */
#define BLAST_SCAN_RANGE_NUCL 0   /**< default scan range (blastn) */

/** length of word to trigger an extension. */
#define BLAST_WORDSIZE_PROT 3   /**< default word size (all protein searches) */
#define BLAST_WORDSIZE_NUCL 11   /**< default word size (blastn) */
#define BLAST_WORDSIZE_MEGABLAST 28   /**< default word size (contiguous 
                                          megablast; for discontig megablast
                                          the word size is explicitly 
                                          overridden) */

#define BLAST_WORDSIZE_MAPPER 18   /**< default word size for mapping rna-seq
                                        to a genome */

/** Default matrix name: BLOSUM62 */
#define BLAST_DEFAULT_MATRIX "BLOSUM62"

/** Protein gap costs are the defaults for the BLOSUM62 scoring matrix.
 * More gap costs are listed in BLASTOptionSetGapParams 
 */

/** cost for the existence of a gap.*/
#define BLAST_GAP_OPEN_PROT 11 /**< default gap open penalty (all 
                                    protein searches) */
#define BLAST_GAP_OPEN_NUCL 5 /**< default gap open penalty (blastn) */
#define BLAST_GAP_OPEN_MEGABLAST 0 /**< default gap open penalty (megablast
                                        with greedy gapped alignment) */
#define BLAST_GAP_OPEN_MAPPER 0

/** cost to extend a gap. */
#define BLAST_GAP_EXTN_PROT 1 /**< default gap open penalty (all 
                                   protein searches) */
#define BLAST_GAP_EXTN_NUCL 2 /**< default gap open penalty (blastn) */
#define BLAST_GAP_EXTN_MEGABLAST 0 /**< default gap open penalty (megablast)
                                        with greedy gapped alignment) */

#define BLAST_GAP_EXTN_MAPPER 4

/** neighboring word score thresholds; a threshold of zero
 *  means that only query and subject words that match exactly
 *  will go into the BLAST lookup table when it is generated 
 */
#define BLAST_WORD_THRESHOLD_BLASTP 11 /**< default neighboring threshold
                                         (blastp and for rpsblast at RPS-BLAST
                                         database creation time) */
#define BLAST_WORD_THRESHOLD_BLASTN 0 /**< default threshold (blastn) */
#define BLAST_WORD_THRESHOLD_BLASTX 12 /**< default threshold (blastx) */
#define BLAST_WORD_THRESHOLD_TBLASTN 13 /**< default neighboring threshold 
                                          (tblastn/rpstblastn) */
#define BLAST_WORD_THRESHOLD_TBLASTX 13 /**< default threshold (tblastx) */
#define BLAST_WORD_THRESHOLD_MEGABLAST 0 /**< default threshold (megablast) */

/** default dropoff for ungapped extension; ungapped extensions
 *  will stop when the score for the extension has dropped from
 *  the current best score by at least this much 
 */
#define BLAST_UNGAPPED_X_DROPOFF_PROT 7 /**< ungapped dropoff score for all
                                             searches except blastn */
#define BLAST_UNGAPPED_X_DROPOFF_NUCL 20 /**< ungapped dropoff score for 
                                              blastn (and megablast) */

/** default dropoff for preliminary gapped extensions */
#define BLAST_GAP_X_DROPOFF_PROT 15 /**< default dropoff (all protein-
                                         based gapped extensions) */
#define BLAST_GAP_X_DROPOFF_NUCL 30 /**< default dropoff for non-greedy
                                         nucleotide gapped extensions */
#define BLAST_GAP_X_DROPOFF_GREEDY 25 /**< default dropoff for greedy
                                         nucleotide gapped extensions */
#define BLAST_GAP_X_DROPOFF_TBLASTX 0 /**< default dropoff for tblastx */

/** default bit score that will trigger gapped extension */
#define BLAST_GAP_TRIGGER_PROT 22.0 /**< default bit score that will trigger
                                         a gapped extension for all protein-
                                         based searches */
#define BLAST_GAP_TRIGGER_NUCL 27.0  /**< default bit score that will trigger
                                         a gapped extension for blastn */

/** default dropoff for the final gapped extension with traceback */
#define BLAST_GAP_X_DROPOFF_FINAL_PROT 25 /**< default dropoff (all protein-
                                               based gapped extensions) */
#define BLAST_GAP_X_DROPOFF_FINAL_NUCL 100 /**< default dropoff for nucleotide
                                               gapped extensions) */
#define BLAST_GAP_X_DROPOFF_FINAL_TBLASTX 0 /**< default dropoff for tblastx */

/** default reward and penalty (only applies to blastn/megablast) */
#define BLAST_PENALTY -3        /**< default nucleotide mismatch score */
#define BLAST_REWARD 1          /**< default nucleotide match score */

#define BLAST_PENALTY_MAPPER -4
#define BLAST_REWARD_MAPPER 1

/** Default parameters for saving hits */
#define BLAST_EXPECT_VALUE 10.0 /**< by default, alignments whose expect
                                     value exceeds this number are discarded */
#define BLAST_HITLIST_SIZE 500 /**< Number of database sequences to save hits 
                                  for */
/** Defaults for PSI-BLAST and DELTA-BLAST options */
#define PSI_INCLUSION_ETHRESH 0.002 /**< Inclusion threshold for PSI BLAST */
#define PSI_PSEUDO_COUNT_CONST 0 /**< Pseudo-count constant for PSI-BLAST */
#define DELTA_INCLUSION_ETHRESH 0.05 /**< Inclusion threshold for DELTA-BLAST */

/** Default genetic code for query and/or database */
#define BLAST_GENETIC_CODE 1  /**< Use the standard genetic code for converting
                                   groups of three nucleotide bases to protein
                                   letters */

/** Default max frequency for a database word. Words with higher frequency
    will be masked in the lookup table. */
#define MAX_DB_WORD_COUNT_MAPPER 30

/** Default maximum insert size: distance on the subject between reads that
    belong to a pair, for spliced and non-spliced alignments */
#define MAGICBLAST_MAX_INSERT_SIZE_SPLICED 1000000
#define MAGICBLAST_MAX_INSERT_SIZE_NONSPLICED 100000


/** Value used to indicate that no IMPALA-style scaling should be performed
 * when scaling a PSSM */
extern const double kPSSM_NoImpalaScaling;

/** Types of the lookup table */
typedef enum {
    eMBLookupTable,  /**< megablast lookup table (includes both
                                contiguous and discontiguous megablast) */
    eSmallNaLookupTable,  /**< lookup table for blastn with small query*/
    eNaLookupTable,  /**< blastn lookup table */
    eAaLookupTable,  /**< standard protein (blastp) lookup table */
    eCompressedAaLookupTable,  /**< compressed alphabet (blastp) lookup table */
    ePhiLookupTable,  /**< protein lookup table specialized for phi-blast */
    ePhiNaLookupTable,  /**< nucleotide lookup table for phi-blast */
    eRPSLookupTable, /**< RPS lookup table (rpsblast and rpstblastn) */
    eIndexedMBLookupTable, /**< use database index as a lookup structure */
    eMixedMBLookupTable, /**< use when some volumes are searched with index and 
                             some are not */
    eNaHashLookupTable  /**< used for 16-base words */
} ELookupTableType;

/** Options needed to construct a lookup table 
 * Also needed: query sequence and query length.
 */
typedef struct LookupTableOptions {
   double threshold; /**< Score threshold for putting words in a lookup table
                          (fractional values are allowed, and could be
                          important if there is scaling involved) */
   ELookupTableType lut_type; /**< What kind of lookup table to construct? */
   Int4 word_size; /**< Determines the size of the lookup table */
   Int4 mb_template_length; /**< Length of the discontiguous words */
   Int4 mb_template_type; /**< Type of a discontiguous word template */
   char* phi_pattern;  /**< PHI-BLAST pattern */
   EBlastProgramType program_number; /**< indicates blastn, blastp, etc. */
   Uint4 stride; /**< number of words to skip after collecting each word */
   bool db_filter; /**< scan the database and include only words that appear
                           in the database between 1 and 9 times
                           (currently implemented only for MB lookuptable
                           and lookup table word size 16) */
   Uint1 max_db_word_count; /**< words with larger frequency in the database
                               will be masked in the lookup table, if the
                               db_filter optoion is on */
} LookupTableOptions;

/** Options for dust algorithm, applies only to nucl.-nucl. comparisons.
 *  value of less than zero means default value will be applied.
 */
typedef struct SDustOptions {
    int level;
    int window;
    int linker;  /**< min distance to link segments. */
} SDustOptions;


/** Options for SEG algorithm, applies only to protein-protein comparisons.
 *  value of less than zero means default value will be applied.
 */
typedef struct SSegOptions {
    int window;     /**< initial window to trigger further work. */
    double locut;
    double hicut;
} SSegOptions;

/// Default value for repeats database filtering
#define kDefaultRepeatFilterDb "repeat/repeat_9606"

/** Filtering options for organsim specific repeats filtering.   
    Currently this consist of only the db name but could be expanded
    in the future to include other types of filtering or other options.
 */
typedef struct SRepeatFilterOptions {
    char* database;   /**< Nucleotide database for mini BLAST search. */
} SRepeatFilterOptions;

/** Filtering options for organism-specific filtering with Window
    Masker.  The taxid and filename are alternative means of choosing
    which Window Masker database to use.
 */
typedef struct SWindowMaskerOptions {
    int          taxid;    /**< Select masking database for this TaxID. */
    const char * database; /**< Use winmasker database at this location. */
} SWindowMaskerOptions;

/** Filtering options for mapping next-generation sequences */
typedef struct SReadQualityOptions {
    double frac_ambig;   /**< Fraction of ambiguous bases */
    int entropy;         /**< Dimer entropy */
} SReadQualityOptions;

/** All filtering options */
typedef struct SBlastFilterOptions {
    bool mask_at_hash;         /**< mask query only for lookup table creation */
    SDustOptions* dustOptions;    /**< low-complexity filtering for nucleotides. */
    SSegOptions* segOptions;      /**< low-complexity filtering for proteins sequences 
            (includes translated nucleotides). */
    SRepeatFilterOptions* repeatFilterOptions;  /**< for organism specific repeat filtering. */
    SWindowMaskerOptions* windowMaskerOptions;  /**< organism specific filtering with window masker. */

    SReadQualityOptions* readQualityOptions;    /**< quality filtering for mapping next-generation sequences */
} SBlastFilterOptions;


/** Options required for setting up the query sequence */
typedef struct QuerySetUpOptions {
   SBlastFilterOptions* filtering_options;  /**< structured options for all filtering 
                             offered from algo/blast/core for BLAST. */
   char* filter_string; /**< DEPRECATED, filtering options above. */

   Uint1 strand_option; /**< In blastn: which strand to search: 1 = forward;
                           2 = reverse; 3 = both */
   Int4 genetic_code;     /**< Genetic code to use for translation, 
                             [t]blastx only */
} QuerySetUpOptions;

/** Options needed for initial word finding and processing */
typedef struct BlastInitialWordOptions {
   double gap_trigger; /**< Score in bits for starting gapped extension */
   Int4 window_size; /**< Maximal allowed distance between 2 hits in case 2 
                        hits are required to trigger the extension */
   Int4 scan_range;  /**< Maximal number of gaps allowed between 2 hits */
   double x_dropoff; /**< X-dropoff value (in bits) for the ungapped 
                         extension */
   EBlastProgramType program_number; /**< indicates blastn, blastp, etc. */
} BlastInitialWordOptions;

/** The algorithm to be used for preliminary
 *  gapped extensions
 */
typedef enum EBlastPrelimGapExt {
    eDynProgScoreOnly,          /**< standard affine gapping */
    eGreedyScoreOnly,           /**< Greedy extension (megaBlast) */
    eJumperWithTraceback,       /**< Jumper extension (mapping) */
    eSmithWatermanScoreOnly     /**< Score-only smith-waterman */
} EBlastPrelimGapExt;

/** The algorithm to be used for final gapped
 *  extensions with traceback
 */
typedef enum EBlastTbackExt {
    eDynProgTbck,          /**< standard affine gapping */
    eGreedyTbck,           /**< Greedy extension (megaBlast) */
    eSmithWatermanTbck,    /**< Smith-waterman finds optimal scores, then 
                                ALIGN_EX to find alignment. */
    eSmithWatermanTbckFull /**< Smith-waterman to find all alignments */
} EBlastTbackExt;

/** Options used for gapped extension 
 *  These include:
 *  a. Penalties for various types of gapping;
 *  b. Drop-off values for the extension algorithms tree exploration;
 *  c. Parameters identifying what kind of extension algorithm(s) should 
 *     be used.
 */
typedef struct BlastExtensionOptions {
   double gap_x_dropoff; /**< X-dropoff value for gapped extension (in bits) */
   double gap_x_dropoff_final;/**< X-dropoff value for the final gapped 
                                  extension (in bits) */
   EBlastPrelimGapExt ePrelimGapExt; /**< type of preliminary gapped extension (normally) for calculating
                              score. */
   EBlastTbackExt eTbackExt; /**< type of traceback extension. */
   Int4 compositionBasedStats; /**< mode of compositional adjustment to use;
                                   if zero then compositional adjustment is
                                   not used */
   Int4 unifiedP; /**< Indicates unified P values to be used in blastp or tblastn */

    Int4 max_mismatches;    /**< Maximum number of mismatches allowed for Jumper */

    Int4 mismatch_window;   /**< Widnow for counting mismatches for Jumper */

   EBlastProgramType program_number; /**< indicates blastn, blastp, etc. */
} BlastExtensionOptions;

/** Options for the Best Hit HSP collection algorithm */
typedef struct BlastHSPBestHitOptions {
    double overhang;
    double score_edge;
} BlastHSPBestHitOptions;

/** Options for the HSP culling algorithm */
typedef struct BlastHSPCullingOptions {
    int max_hits; /**< Maximum number of hits per area of query. */
} BlastHSPCullingOptions;

typedef struct BlastHSPSubjectBestHitOptions {
    unsigned int max_range_diff;
} BlastHSPSubjectBestHitOptions;

/** Structure containing the HSP filtering/writing options */
typedef struct BlastHSPFilteringOptions {
    /** Best Hit algorithm */
    BlastHSPBestHitOptions* best_hit;
    EBlastStage best_hit_stage; /*<< when to apply the best hit algorithm */

    /** culling algorithm */
    BlastHSPCullingOptions* culling_opts;
    EBlastStage culling_stage; /*<< when to apply the culling algorithm */

    /** Subject Culling */
    BlastHSPSubjectBestHitOptions * subject_besthit_opts;
} BlastHSPFilteringOptions;

/** Options used when evaluating and saving hits
 *  These include: 
 *  a. Restrictions on the number of hits to be saved;
 *  b. Restrictions on the quality and positions of hits to be saved;
 *  c. Parameters used to evaluate the quality of hits.
 */
typedef struct BlastHitSavingOptions {
   double expect_value; /**< The expect value cut-off threshold for an HSP, or
                            a combined hit if sum statistics is used */
   Int4 cutoff_score; /**< The (raw) score cut-off threshold */
   Int4 cutoff_score_fun[2]; /**< Coefficients x100 for the raw score cut-off
                                 threshold as a function of query length:
                                 x[0] + x[1] * query_length*/
   double percent_identity; /**< The percent identity cut-off threshold */

   Int4 max_edit_distance; /**< Maximum number of mismatches and gaps */

   Int4 hitlist_size;/**< Maximal number of database sequences to return
                        results for */
   Int4 hsp_num_max; /**< Maximal number of HSPs to save for one database 
                        sequence */
   Int4 total_hsp_limit; /**< Maximal total number of HSPs to keep */
   Int4 culling_limit; /**< If the query range of an HSP is contained in
                            at least this many higher-scoring HSPs, throw
                            away the HSP as redundant (turned off if zero) */
   Int4 mask_level; /**< Only keep the highest scoring HSP when more than
                          one HSP overlaps the same region of the query by
                          more than or equal to mask_level %. -RMH- */

   /********************************************************************/
   /* Merge all these in a structure for clarity? */
   bool do_sum_stats; /**< Force sum statistics to be used to combine HSPs,
                          TRUE by default for all ungapped searches and translated
                          gapped searches (except RPS-BLAST) */
   Int4 longest_intron; /**< The longest distance between HSPs allowed for
                           combining via sum statistics with uneven gaps */
   /********************************************************************/

   Int4 min_hit_length;    /**< optional minimum alignment length; alignments
                                not at least this long are discarded */
   Int4 min_diag_separation; /**< How many diagonals separate a hit from a substantial alignment
                                  before it's not blocked out. Must be > 0 to be used. */
   EBlastProgramType program_number; /**< indicates blastn, blastp, etc. */

   /** Contains options to configure the HSP filtering/writering structures
    * If not set, the default HSP filtering strategy is used.
    */
   BlastHSPFilteringOptions* hsp_filt_opt;

   /** Low-score option.  Do not pass ungapped alignments on for later processing if
    * the hitlist is already full of other alignments unless the ungapped aligment 
    * is above the fraction X of the least significant database match.
    * zero should turn this off.
    */
   double low_score_perc;

   double query_cov_hsp_perc; /**< Min query coverage hsp percentage */

    /* Used by default hsp filtering strategy, num of best hsps to keep per subject
     * seq for each query. Note that hsp_num_max should be used only to reduce memory footprint,
     * it does not guarantee best hsp per query due to query concatenation
     */
   Int4 max_hsps_per_subject;

   /**< Queries are paired reads, for mapping */
   bool paired;
   /**< Splice HSPs for each query (for mapping RNA-Seq to a genome) */
   bool splice;

} BlastHitSavingOptions;

/** Scoring options block 
 *  Used to produce the BlastScoreBlk structure
 *  This structure may be needed for lookup table construction (proteins only),
 *  and for evaluating alignments. 
 */
typedef struct BlastScoringOptions {
   char* matrix;   /**< Name of the matrix containing all scores: needed for
                        finding neighboring words */
   char* matrix_path; /**< Directory path to where matrices are stored. */
   Int2 reward;      /**< Reward for a match */
   Int2 penalty;     /**< Penalty for a mismatch */
   bool gapped_calculation; /**< gap-free search if FALSE */
   bool complexity_adjusted_scoring; /**< Use cross_match-like complexity
                                           adjustment on raw scores. -RMH- */
   Int4 gap_open;    /**< Extra penalty for starting a gap */
   Int4 gap_extend;  /**< Penalty for each gap residue */

   /* only blastx and tblastn (When query & subj are diff) */
   bool is_ooframe; /**< Should out-of-frame gapping be used in a translated
                          search? */
   Int4 shift_pen;   /**< Penalty for shifting a frame in out-of-frame 
                        gapping */
   EBlastProgramType program_number; /**< indicates blastn, blastp, etc. */
} BlastScoringOptions;

/** Options for setting up effective lengths and search spaces.  
 * The values are those the user has specified to override the real sizes.
 */
typedef struct BlastEffectiveLengthsOptions {
   Int8 db_length;    /**< Database length to be used for statistical
                         calculations */
   Int4 dbseq_num;    /**< Number of database sequences to be used for
                           statistical calculations */
   Int4 num_searchspaces; /**< Number of elements in searchsp_eff, this must be
                            equal to the number of contexts in the search */
   Int8 *searchsp_eff; /**< Search space to be used for statistical
                           calculations (one such per query context) */
} BlastEffectiveLengthsOptions;

/** Options used in protein BLAST only (PSI, PHI, RPS and translated BLAST)
 *  Some of these possibly should be transfered elsewhere  
 */
typedef struct PSIBlastOptions {
    /** Pseudocount constant. Needed for the computing the PSSM residue 
     * frequencies */
    Int4 pseudo_count;

    /*** The following options are used at the API layer to specify how the
     * multiple sequence alignment is built from pairwise alignments. These
     * could go in their own structure in the future. */

    /** Minimum evalue for inclusion in PSSM calculation. Needed for the
     * conversion of Seq-align into a multiple sequence alignment and for
     * composition based statistics */
    double inclusion_ethresh;

    /** If set to TRUE, use the best alignment when multiple HSPs are found 
     * in a query-subject alignment (i.e.: HSP with the lowest e-value), else
     * use all HSPs in a query-subject alignment. This option does not apply to
     * the PSSM engine, it applies to the processing of pairwise sequence
     * alignments to build a multiple sequence alignment structure 
     * @sa CPsiBlastInputData (to be implemented)
     */
    bool use_best_alignment;

    /** Compatibility option for the NCBI's structure group (note 
     * nsg_ prefix, stands for NCBI's structure group). When set to true, the
     * PSSM engine will function in the same way the C toolkit PSSM engine did
     * in the structure group's cddumper application. This option should be 
     * set to FALSE by default as it enables the following behavior in the 
     * PSSM engine:
     * <pre>
     * 1) Ignores the query sequence (on certain stages of PSSM creation only)
     * 2) Skips validation of multiple sequence alignment data
     * 3) Disables assertions and validation in _PSICheckSequenceWeights
     * 4) If no aligned sequences are provided in the multiple sequence
     * alignment, NULL PSSM frequency ratios are returned and the PSSM is built
     * based on the underlying substitution matrix.
     * </pre>
     * Do not set this to TRUE unless you know what you are doing.
     */
    bool nsg_compatibility_mode;

    /** Scaling factor as used in IMPALA to do the matrix rescaling. Default
     * value of 1.0 means not to use it. Makemat/formatrpsdb set this value to
     * 100 by default, Kappa_RedoAlignmentCore uses 32. Provided so that the
     * NCBI structure group can create scaled PSSMs as the output of the PSSM
     * engine. Do not change this unless you know what you are doing.
     */
    double impala_scaling_factor;

    /** This turns off a validation for the multiple sequence alignment in the
     * PSSM engine for unaligned positions. Needed when a multiple sequence
     * alignment is provided on the command line (e.g.: -in_msa option). 
     */
    bool ignore_unaligned_positions;

} PSIBlastOptions;


/** Options used to create the ReadDBFILE structure 
 *  Include database name and various information for restricting the database
 *  to a subset.
 */
typedef struct BlastDatabaseOptions {
   Int4 genetic_code;   /**< Genetic code to use for translation, 
                             tblast[nx] only */
} BlastDatabaseOptions;

/********************************************************************************

    Functions to create options blocks with default values
    and free them after use.

*********************************************************************************/

/** Frees SDustOptions.
 * @param dust_options object to free
 * @return NULL pointer
 */
NCBI_XBLAST_EXPORT
SDustOptions* SDustOptionsFree(SDustOptions* dust_options);

/** Allocates memory for SDustOptions, fills in defaults.
 * @param dust_options options that are being returned [in|out]
 * @return zero on sucess
 */
NCBI_XBLAST_EXPORT
Int2 SDustOptionsNew(SDustOptions* *dust_options);

/** Frees SSegOptions.
 * @param seg_options object to free [in]
 * @return NULL pointer
 */
NCBI_XBLAST_EXPORT
SSegOptions* SSegOptionsFree(SSegOptions* seg_options);

/** Allocates memory for SSegOptions, fills in defaults. [in|out]
 * @param seg_options options that are being returned [in|out]
 * @return zero on sucess
 */
NCBI_XBLAST_EXPORT
Int2 SSegOptionsNew(SSegOptions* *seg_options);

/** Resets name of db for repeat filtering.
 * @param repeat_options already allocated options constaining field to be reset [in|out]
 * @param dbname name of the database(s) [in]
 * @return zero on sucess
 */
NCBI_XBLAST_EXPORT
Int2 SRepeatFilterOptionsResetDB(SRepeatFilterOptions* *repeat_options, const char* dbname);

/** Resets name of db for window masker filtering.
 * @param winmask_options options block constaining field to be reset [in|out]
 * @param dbname name of the database(s) [in]
 * @return zero on sucess
 */
NCBI_XBLAST_EXPORT
Int2 SWindowMaskerOptionsResetDB(SWindowMaskerOptions ** winmask_options,
                                 const char            * dbname);

/** Frees SRepeatFilterOptions.
 * @param repeat_options object to free [in]
 * @return NULL pointer
 */
NCBI_XBLAST_EXPORT
SRepeatFilterOptions* SRepeatFilterOptionsFree(SRepeatFilterOptions* repeat_options);

/** Frees SWindowMaskerOptions.
 * @param winmask_options object to free [in]
 * @return NULL pointer
 */
NCBI_XBLAST_EXPORT SWindowMaskerOptions*
SWindowMaskerOptionsFree(SWindowMaskerOptions * winmask_options);

/** Allocates memory for SRepeatFilterOptions, fills in defaults.
 * @param repeat_options options that are being returned [in|out]
 * @return zero on sucess
 */
NCBI_XBLAST_EXPORT
Int2 SRepeatFilterOptionsNew(SRepeatFilterOptions ** repeat_options);

/** Allocates memory for SWindowMaskerOptions, fills in defaults.
 * @param winmask_options options that are being returned [in|out]
 * @return zero on sucess
 */
NCBI_XBLAST_EXPORT
Int2 SWindowMaskerOptionsNew(SWindowMaskerOptions ** winmask_options);

/** Allocates memory for SReadQualityOptions, fills in defaults.
 * @param read_quality_ptions options that are being returned [in|out]
 * @return zero on sucess
 */
NCBI_XBLAST_EXPORT
Int2 SReadQualityOptionsNew(SReadQualityOptions ** read_quality_options);

/** Frees memory for SReadQualityOptions */
NCBI_XBLAST_EXPORT
SReadQualityOptions* SReadQualityOptionsFree(
                                 SReadQualityOptions * read_quality_options);


/** Frees SBlastFilterOptions and all subservient structures.
 * @param filter_options object to free
 * @return NULL pointer
 */
NCBI_XBLAST_EXPORT
SBlastFilterOptions* SBlastFilterOptionsFree(SBlastFilterOptions* filter_options);

/**  Merges two sets of options together, taking the non-default one as preferred.  if
 * both are non-default then one or the other is taken.
 * @param combined object that is returned [in|out]
 * @param opt1 first set of options [in]
 * @param opt2 second set of options [in]
 * @return zero on success. 
 */
NCBI_XBLAST_EXPORT
Int2 SBlastFilterOptionsMerge(SBlastFilterOptions** combined, const SBlastFilterOptions* opt1,
       const SBlastFilterOptions* opt2);

/** Types of filtering options. */
typedef enum EFilterOptions {
    eSeg,            /**< low-complexity for proteins. */
    eDust,           /**< low-complexity for nucleotides. */
    eRepeats,         /**< Repeat filtering for nucleotides. */
    eDustRepeats,    /**< Repeat and dust filtering for nucleotides. */
    eEmpty       /**< no filtering at all. */
} EFilterOptions;

/** Allocates memory for SBlastFilterOptions and 
 * @param filter_options options that are being returned [in|out]
 * @param type specify either dust or seg (now) with EFilterOptions [in]
 * @return zero on sucess
 */
NCBI_XBLAST_EXPORT
Int2 SBlastFilterOptionsNew(SBlastFilterOptions* *filter_options, EFilterOptions type);

/** Queries whether no masking is required
 * @param filter_options the object to be queried [in]
 * @return TRUE if no filtering is required or argument is NULL, FALSE
 * otherwise
 */
NCBI_XBLAST_EXPORT
bool SBlastFilterOptionsNoFiltering(const SBlastFilterOptions* filter_options);

/** Queries whether masking should be done only for the lookup table or for the entire search.
 * @param filter_options the object to be queried [in]
 * @return TRUE or FALSE, FALSE if filter_options is NULL.
 */
NCBI_XBLAST_EXPORT
bool SBlastFilterOptionsMaskAtHash(const SBlastFilterOptions* filter_options);

/** Validates filter options to ensure that program and options are consistent
 * and that options have valid values.
 * @param program_number Program number (blastn, blastp, etc.) [in]
 * @param filter_options options to add to [in]
 * @param blast_message error or warning (optional) [out] 
 * @return zero on success
 */
NCBI_XBLAST_EXPORT
Int2 SBlastFilterOptionsValidate(EBlastProgramType program_number, const SBlastFilterOptions* filter_options, 
       Blast_Message* *blast_message);


/** Deallocate memory for QuerySetUpOptions. 
 * @param options Structure to free [in]
 */
NCBI_XBLAST_EXPORT
QuerySetUpOptions* BlastQuerySetUpOptionsFree(QuerySetUpOptions* options);


/** Allocate memory for QuerySetUpOptions and fill with default values.  
 * @param options The options that have are being returned [out]
 */
NCBI_XBLAST_EXPORT
Int2 BlastQuerySetUpOptionsNew(QuerySetUpOptions* *options);

/** Fill non-default contents of the QuerySetUpOptions.
 * @param options The options structure [in] [out]  
 * @param program Program number (blastn, blastp, etc.) [in]
 * @param filter_string Parsable string of filtering options [in] 
 * @param strand_option which strand to search [in]
*/
NCBI_XBLAST_EXPORT
Int2 BLAST_FillQuerySetUpOptions(QuerySetUpOptions* options,
        EBlastProgramType program, const char *filter_string, Uint1 strand_option);


/** Deallocate memory for BlastInitialWordOptions.
 * @param options Structure to free [in]
 */
NCBI_XBLAST_EXPORT
BlastInitialWordOptions*
BlastInitialWordOptionsFree(BlastInitialWordOptions* options);

/** Allocate memory for BlastInitialWordOptions and fill with default values.
 * @param program Program number (blastn, blastp, etc.) [in]
 * @param options The options that have are being returned [out] 
*/
NCBI_XBLAST_EXPORT
Int2
BlastInitialWordOptionsNew(EBlastProgramType program, 
   BlastInitialWordOptions* *options);

/** Validate correctness of the initial word options.
 * @param program_number Type of BLAST program [in]
 * @param options Initial word options [in]
 * @param blast_msg Describes any validation problems found [out]
 * @return Validation status
 */
NCBI_XBLAST_EXPORT
Int2
BlastInitialWordOptionsValidate(EBlastProgramType program_number,
   const BlastInitialWordOptions* options, 
   Blast_Message* *blast_msg);

/** Fill non-default values in the BlastInitialWordOptions structure.
 * @param options The options structure [in] [out] 
 * @param program Program number (blastn, blastp, etc.) [in]
 * @param window_size Size of a largest window between 2 words for the two-hit
 *                    version [in]
 * @param xdrop_ungapped The value of the X-dropoff for ungapped extensions [in]
*/
NCBI_XBLAST_EXPORT
Int2
BLAST_FillInitialWordOptions(BlastInitialWordOptions* options, 
                EBlastProgramType program, 
                Int4 window_size, double xdrop_ungapped);

/** Deallocate memory for BlastExtensionOptions.
 * @param options Structure to free [in]
 */
NCBI_XBLAST_EXPORT
BlastExtensionOptions*
BlastExtensionOptionsFree(BlastExtensionOptions* options);

/** Allocate memory for BlastExtensionOptions and fill with default values.
 * @param program Program number (blastn, blastp, etc.) [in]
 * @param options The options that are being returned [out]
 * @param gapped The search is gapped [in]
*/
NCBI_XBLAST_EXPORT
Int2
BlastExtensionOptionsNew(EBlastProgramType program, BlastExtensionOptions* *options, bool gapped);

/** Fill non-default values in the BlastExtensionOptions structure.
 * @param options The options structure [in] [out]
 * @param program Program number (blastn, blastp, etc.) [in]
 * @param greedy In how many stages of the search greedy alignment is 
 *               used (values 0, 1, 2)? FIXME [in]
 * @param x_dropoff X-dropoff parameter value for preliminary gapped 
 *                  extensions [in]
 * @param x_dropoff_final X-dropoff parameter value for final gapped 
 *                        extensions with traceback [in]
 * @todo the greedy parameter to this function is tied to the blast_driver's
 * command line argument for greedy... couldn't this be EBlastPrelimGapExt?
*/
NCBI_XBLAST_EXPORT
Int2
BLAST_FillExtensionOptions(BlastExtensionOptions* options, 
   EBlastProgramType program, Int4 greedy, double x_dropoff, 
   double x_dropoff_final);



/** Validate contents of BlastExtensionOptions.
 * @param program_number Type of BLAST program [in]
 * @param options Options to be validated [in]
 * @param blast_msg Describes any validation problems found [out]
*/
NCBI_XBLAST_EXPORT
Int2 BlastExtensionOptionsValidate(EBlastProgramType program_number, 
        const BlastExtensionOptions* options, Blast_Message* *blast_msg);

/**  Deallocate memory for BlastScoringOptions. 
 * @param options Structure to free [in]
 */
NCBI_XBLAST_EXPORT
BlastScoringOptions* BlastScoringOptionsFree(BlastScoringOptions* options);

/** Allocate memory for BlastScoringOptions and fill with default values. 
 * @param program Program number (blastn, blastp, etc.) [in]
 * @param options The options that are being returned [out]
*/
NCBI_XBLAST_EXPORT
Int2 BlastScoringOptionsNew(EBlastProgramType program, BlastScoringOptions* *options);

/** Fill non-default values in the BlastScoringOptions structure. 
 * @param options The options structure [in] [out]
 * @param program Program number (blastn, blastp, etc.) [in]
 * @param greedy_extension Is greedy extension algorithm used? [in]
 * @param penalty Mismatch penalty score (blastn only) [in]
 * @param reward Match reward score (blastn only) [in]
 * @param matrix Name of the BLAST matrix (all except blastn) [in]
 * @param gap_open Extra cost for opening a gap [in]
 * @param gap_extend Cost of a gap [in]
*/
NCBI_XBLAST_EXPORT
Int2 
BLAST_FillScoringOptions(BlastScoringOptions* options, EBlastProgramType program, 
   bool greedy_extension, Int4 penalty, Int4 reward, const char *matrix, 
   Int4 gap_open, Int4 gap_extend);

/** Validate contents of BlastScoringOptions.
 * @param program_number Type of BLAST program [in]
 * @param options Options to be validated [in]
 * @param blast_msg Describes any validation problems found [out]
*/
NCBI_XBLAST_EXPORT
Int2
BlastScoringOptionsValidate(EBlastProgramType program_number, 
   const BlastScoringOptions* options, Blast_Message* *blast_msg);

/** Produces copy of "old" options, with new memory allocated.
 * @param new_opt Contains copied BlastScoringOptions upon return [out]
 * @param old_opt BlastScoringOptions to be copied [in]
*/
NCBI_XBLAST_EXPORT
Int2 BlastScoringOptionsDup(BlastScoringOptions* *new_opt, const BlastScoringOptions* old_opt);

/** Resets matrix name option. Automatically converts the name to upper case.
 * @param opts Options structure to update. [in] [out]
 * @param matrix_name New matrix name. If NULL, old matrix name is left 
 *                    as is. [in]
 */
NCBI_XBLAST_EXPORT
Int2 BlastScoringOptionsSetMatrix(BlastScoringOptions* opts,
                                  const char* matrix_name);


/** Deallocate memory for BlastEffectiveLengthsOptions*. 
 * @param options Structure to free [in]
 */
NCBI_XBLAST_EXPORT
BlastEffectiveLengthsOptions* 
BlastEffectiveLengthsOptionsFree(BlastEffectiveLengthsOptions* options);

/** Allocate memory for BlastEffectiveLengthsOptions* and fill with 
 * default values. 
 * @param options The options that are being returned [out]
 */
NCBI_XBLAST_EXPORT
Int2 BlastEffectiveLengthsOptionsNew(BlastEffectiveLengthsOptions* *options);

/** Return true if the search spaces is set for any of the queries in the
 * search
 * @param options The options to examine [in]
 */
NCBI_XBLAST_EXPORT
bool
BlastEffectiveLengthsOptions_IsSearchSpaceSet(const
                                              BlastEffectiveLengthsOptions*
                                              options);

/** Fill the non-default values in the BlastEffectiveLengthsOptions structure.
 * @param options The options [in] [out]
 * @param dbseq_num Number of sequences in the database (if zero real value will be used) [in]
 * @param db_length Total length of the database (if zero real value will be used) [in]
 * @param *searchsp_eff Array of effective search spaces (the real value 
 *                   will be used for elements that are 0). If array 
 *                   contains one element, all contexts use this value. 
 *                   If array has multiple elements, the number must match
 *                   the number of contexts in the search [in]
 * @param num_searchsp The number of elements in searchsp_eff [in]
 */
NCBI_XBLAST_EXPORT
Int2 
BLAST_FillEffectiveLengthsOptions(BlastEffectiveLengthsOptions* options, 
                                  Int4 dbseq_num, Int8 db_length, 
                                  Int8 *searchsp_eff, Int4 num_searchsp);


/** Allocate memory for lookup table options and fill with default values.
 * @param program Program number (blastn, blastp, etc.) [in]
 * @param options The options that are being returned [out]
 */
NCBI_XBLAST_EXPORT
Int2 LookupTableOptionsNew(EBlastProgramType program, LookupTableOptions* *options);


/** Allocate memory for lookup table options and fill with default values.
 * @param options The options [in] [out]
 * @param program Program number (blastn, blastp, etc.) [in]
 * @param is_megablast Megablast (instead of blastn) if TRUE [in]
 * @param threshold Threshold value for finding neighboring words
                    (fractional values are allowed, though unless
                    the engine scales up alignment scores a fractional
                    threshold will be rounded down) [in]
 * @param word_size Number of matched residues in an initial word [in]
 */
NCBI_XBLAST_EXPORT
Int2 
BLAST_FillLookupTableOptions(LookupTableOptions* options, 
   EBlastProgramType program, bool is_megablast, double threshold,
   Int4 word_size);


/** Deallocates memory for LookupTableOptions*.
 * @param options Structure to free [in]
 */
NCBI_XBLAST_EXPORT
LookupTableOptions*
LookupTableOptionsFree(LookupTableOptions* options);

/** Validate LookupTableOptions.
 * @param program_number BLAST program [in]
 * @param options The options that have are being returned [in]
 * @param blast_msg Describes any validation problems found [out]
*/
NCBI_XBLAST_EXPORT
Int2
LookupTableOptionsValidate(EBlastProgramType program_number, 
   const LookupTableOptions* options,  Blast_Message* *blast_msg);

/** Deallocate memory for BlastHitSavingOptions. 
 * @param options Structure to free [in]
 */
NCBI_XBLAST_EXPORT
BlastHitSavingOptions*
BlastHitSavingOptionsFree(BlastHitSavingOptions* options);

/** Validate BlastHitSavingOptions
 * @param program_number BLAST program [in]
 * @param options The options that have are being returned [in]
 * @param blast_msg Describes any validation problems found [out]
*/

NCBI_XBLAST_EXPORT
Int2
BlastHitSavingOptionsValidate(EBlastProgramType program_number,
   const BlastHitSavingOptions* options, Blast_Message* *blast_msg);

/** Allocate memory for BlastHitSavingOptions.
 * @param program Program number (blastn, blastp, etc.) [in]
 * @param options The options that are being returned [out]
 * @param gapped_calculation is this search gapped? [in]
*/
NCBI_XBLAST_EXPORT
Int2 BlastHitSavingOptionsNew(EBlastProgramType program, 
        BlastHitSavingOptions** options,
        bool gapped_calculation);

/** Allocate memory for BlastHitSavingOptions.
 * @param options The options [in] [out]
 * @param evalue The expected value threshold [in]
 * @param hitlist_size How many database sequences to save per query? [in]
 * @param is_gapped is this a gapped alignment? [in]
 * @param culling_limit Number of higher-scoring HSPs that must contain
 *                      the query range of an HSP before that HSP is declared
 *                      to be redundant (ignored if zero) [in]
 * @param min_diag_separation Delete HSPs whose endpoints are at most this 
 *                       many diagonals from a higher-scoring HSP. If zero,
 *                       delete HSPs whose query and subject ranges are
 *                       enveloped by those of a higher-scoring HSP [in]
*/
NCBI_XBLAST_EXPORT
Int2
BLAST_FillHitSavingOptions(BlastHitSavingOptions* options, 
                           double evalue, Int4 hitlist_size,
                           bool is_gapped,
                           Int4 culling_limit,
                           Int4 min_diag_separation);

/** Initialize default options for PSI BLAST 
 * @param psi_options pointer to pointer where structure will be allocated [in]
 * @return 1 in case of memory allocation failure or if psi_options is NULL, 0
 * in case of success
 */
NCBI_XBLAST_EXPORT
Int2 PSIBlastOptionsNew(PSIBlastOptions** psi_options);

/** Validates the PSI BLAST options so that they have sane values.
 * @param psi_options structure to validate [in]
 * @param blast_msg Describes any validation problems found [out]
 * @return 0 on success 1 on failure
 */
Int2 PSIBlastOptionsValidate(const PSIBlastOptions* psi_options,
                             Blast_Message** blast_msg);

/** Deallocate PSI BLAST options */
NCBI_XBLAST_EXPORT
PSIBlastOptions* PSIBlastOptionsFree(PSIBlastOptions* psi_options);

/** Allocate and initialize a BlastHSPBestHitOptions structure */
NCBI_XBLAST_EXPORT
BlastHSPBestHitOptions* BlastHSPBestHitOptionsNew(double overhang, 
                                                  double score_edge);

/** Validate the best hit algorithm parameters (if any) in the
 * @param opts BlastHSPFilteringOptions structure 
 * @return 0 on success, else non-zero
 */
NCBI_XBLAST_EXPORT
Int2
BlastHSPBestHitOptionsValidate(const BlastHSPFilteringOptions* opts);

/** Deallocate a BlastHSPBestHitOptions structure 
 * @param opt object to be deallocated. [in]
 */
NCBI_XBLAST_EXPORT
BlastHSPBestHitOptions* BlastHSPBestHitOptionsFree(BlastHSPBestHitOptions* opt);

/** Allocate a new object for culling options.
 * @param max number of HSPs that may be aligned to one part of query [in]
 */
NCBI_XBLAST_EXPORT
BlastHSPCullingOptions* BlastHSPCullingOptionsNew(int max);

/** Validate culling options.
 * @param opts BlastHSPFilteringOptions structure 
 * @return 0 on success, else non-zero
 */
NCBI_XBLAST_EXPORT
Int2
BlastHSPCullingOptionsValidate(const BlastHSPFilteringOptions* opts);

/** Deallocates culling options structure.
 * @param culling_opts object to be deallocated. [in]
 */
NCBI_XBLAST_EXPORT
BlastHSPCullingOptions* 
BlastHSPCullingOptionsFree(BlastHSPCullingOptions* culling_opts);

/** Allocate and initialize a BlastHSPFilteringOptions structure */
NCBI_XBLAST_EXPORT
BlastHSPFilteringOptions* BlastHSPFilteringOptionsNew();

/** Add the best hit options. Responsibility for best_hit is taken over by the 
 * BlastHSPFilteringOptions
 * @param filt_opts HSP filtering options [in]
 * @param best_hit Best Hit algorithm options. Ownership of this is taken by
 * the BlastHSPFilteringOptions structure [in|out]
 */
NCBI_XBLAST_EXPORT
Int2
BlastHSPFilteringOptions_AddBestHit(BlastHSPFilteringOptions* filt_opts,
                                    BlastHSPBestHitOptions** opts, 
                                    EBlastStage stage);
/** Validates the BlastHSPFilteringOptions structure */
NCBI_XBLAST_EXPORT
Int2
BlastHSPFilteringOptions_AddCulling(BlastHSPFilteringOptions* filt_opts,
                                    BlastHSPCullingOptions** opts, 
                                    EBlastStage stage);

/** Validates the BlastHSPFilteringOptions structure */
NCBI_XBLAST_EXPORT
Int2
BlastHSPFilteringOptionsValidate(const BlastHSPFilteringOptions* opts);

/** Deallocate a BlastHSPFilteringOptions structure */
NCBI_XBLAST_EXPORT
BlastHSPFilteringOptions*
BlastHSPFilteringOptionsFree(BlastHSPFilteringOptions* opts);

/** Allocates the BlastDatabase options structure and sets the default
 * database genetic code value (BLAST_GENETIC_CODE). Genetic code string in
 * ncbistdaa must be populated by client code */
NCBI_XBLAST_EXPORT
Int2 BlastDatabaseOptionsNew(BlastDatabaseOptions** db_options);

/** Deallocate database options */
NCBI_XBLAST_EXPORT
BlastDatabaseOptions* 
BlastDatabaseOptionsFree(BlastDatabaseOptions* db_options);

/** Initialize all the BLAST search options structures with the default
 * values.
 * @param blast_program Type of blast program: blastn, blastp, blastx, 
 *                      tblastn, tblastx) [in]
 * @param lookup_options Lookup table options [out]
 * @param query_setup_options Query options [out]
 * @param word_options Initial word processing options [out]
 * @param ext_options Extension options [out]
 * @param hit_options Hit saving options [out]
 * @param score_options Scoring options [out]
 * @param eff_len_options Effective length options [out]
 * @param protein_options Protein BLAST options [out]
 * @param db_options BLAST database options [out]
 */
NCBI_XBLAST_EXPORT
Int2 BLAST_InitDefaultOptions(EBlastProgramType blast_program,
   LookupTableOptions** lookup_options,
   QuerySetUpOptions** query_setup_options, 
   BlastInitialWordOptions** word_options,
   BlastExtensionOptions** ext_options,
   BlastHitSavingOptions** hit_options,
   BlastScoringOptions** score_options,
   BlastEffectiveLengthsOptions** eff_len_options,
   PSIBlastOptions** protein_options,
   BlastDatabaseOptions** db_options);

/** Validate all options */
NCBI_XBLAST_EXPORT
Int2 BLAST_ValidateOptions(EBlastProgramType program_number,
                           const BlastExtensionOptions* ext_options,
                           const BlastScoringOptions* score_options, 
                           const LookupTableOptions* lookup_options, 
                           const BlastInitialWordOptions* word_options, 
                           const BlastHitSavingOptions* hit_options,
                           Blast_Message* *blast_msg);



/** Get thresholds for word-finding suggested by Stephen Altschul.
 *
 * @param program_number Type of blast program: blastn, blastp, blastx, 
 *                      tblastn, tblastx) [in]
 * @param matrixName matrix, e.g., BLOSUM62 [in]
 * @param threshold returns suggested value [in|out]
 * @return zero on success
 */
NCBI_XBLAST_EXPORT
Int2 BLAST_GetSuggestedThreshold(EBlastProgramType program_number, 
                                 const char* matrixName, 
                                 double* threshold);

/** Get window sizes for two hit algorithm suggested by Stephen Altschul.
 *
 * @param program_number Type of blast program: blastn, blastp, blastx, 
 *                      tblastn, tblastx) [in]
 * @param matrixName matrix, e.g., BLOSUM62 [in]
 * @param window_size returns suggested value [in|out]
 * @return zero on success
 */
NCBI_XBLAST_EXPORT
Int2 BLAST_GetSuggestedWindowSize(EBlastProgramType program_number, 
                                 const char* matrixName, 
                                 Int4* window_size);

/** Allocate a new object for subject besthit options.
 * @params isProtein true if protein alignment [in]
 */
NCBI_XBLAST_EXPORT
BlastHSPSubjectBestHitOptions* BlastHSPSubjectBestHitOptionsNew(bool isProtein);

/** Validate subject besthit options.
 * @param opts BlastHSPFilteringOptions structure
 * @return 0 on success, else non-zero
 */
NCBI_XBLAST_EXPORT
Int2
BlastHSPSubjectBestHitOptionsValidate(const BlastHSPFilteringOptions* opts);

/** Deallocates subject besthit structure.
 * @param subject_besthit object to be deallocated. [in]
 */
NCBI_XBLAST_EXPORT
BlastHSPSubjectBestHitOptions*
BlastHSPSubjectBestHitOptionsFree(BlastHSPSubjectBestHitOptions* subject_besthit_opts);

NCBI_XBLAST_EXPORT
Int2
BlastHSPFilteringOptions_AddSubjectBestHit(BlastHSPFilteringOptions* filt_opts,
                                           BlastHSPSubjectBestHitOptions** subject_besthit);

#define DEFAULT_SUBJECT_BESTHIT_PROT_MAX_RANGE_DIFF 3
#define DEFAULT_SUBJECT_BESTHIT_NUCL_MAX_RANGE_DIFF 3

#ifdef __cplusplus
}
#endif
#endif /* !__BLASTOPTIONS__ */

