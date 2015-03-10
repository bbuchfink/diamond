/* $Id: blast_options.c 406210 2013-07-11 14:33:14Z madden $
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

/** @file blast_options.c
 *  The structures and functions in blast_options.[ch] should be used to specify 
 *  user preferences.  The options structures should not be changed by the BLAST code
 *  but rather be read to determine user preferences.  When possible these structures
 *  should be passed in as "const".
 *
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: blast_options.c 406210 2013-07-11 14:33:14Z madden $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include "blast_options.h"
#include "blast_filter.h"
#include "blast_stat.h"

const double kPSSM_NoImpalaScaling = 1.0;

/** Declared in blast_def.h as extern const. */
const int kDustLevel = 20;
const int kDustWindow = 64;
const int kDustLinker = 1;

SDustOptions* SDustOptionsFree(SDustOptions* dust_options)
{
    if (dust_options)
      sfree(dust_options);
    return NULL;
}

Int2 SDustOptionsNew(SDustOptions* *dust_options)
{
    if (dust_options == NULL)
        return 1;

    *dust_options = (SDustOptions*) malloc(sizeof(SDustOptions));
    (*dust_options)->level = kDustLevel;
    (*dust_options)->window = kDustWindow;
    (*dust_options)->linker = kDustLinker;

    return 0;
}

SSegOptions* SSegOptionsFree(SSegOptions* seg_options)
{
    if (seg_options)
      sfree(seg_options);
    return NULL;
}

Int2 SSegOptionsNew(SSegOptions* *seg_options)
{
    if (seg_options == NULL)
        return 1;

    *seg_options = (SSegOptions*) malloc(sizeof(SSegOptions));
    (*seg_options)->window = kSegWindow;
    (*seg_options)->locut = kSegLocut;
    (*seg_options)->hicut = kSegHicut;

    return 0;
}

Int2 SWindowMaskerOptionsNew(SWindowMaskerOptions ** winmask_options)
{
    if (winmask_options) {
        *winmask_options = (SWindowMaskerOptions*) calloc(1, sizeof(SWindowMaskerOptions));
        if (*winmask_options == NULL)
            return BLASTERR_MEMORY;
        
        (*winmask_options)->taxid = 0;
        (*winmask_options)->database = NULL;
        return 0;
    }
    return 1;
}

SWindowMaskerOptions* SWindowMaskerOptionsFree(SWindowMaskerOptions* winmask_options)
{
    if (winmask_options)
    {
        if (winmask_options->database)
        {
            sfree(winmask_options->database);
        }
        sfree(winmask_options);
    }
    return NULL;
}

SRepeatFilterOptions* SRepeatFilterOptionsFree(SRepeatFilterOptions* repeat_options)
{
    if (repeat_options)
    {
        sfree(repeat_options->database);
        sfree(repeat_options);
    }
    return NULL;
}

Int2 SRepeatFilterOptionsNew(SRepeatFilterOptions* *repeat_options)
{

    if (repeat_options == NULL)
        return 1;

    *repeat_options = (SRepeatFilterOptions*) calloc(1, sizeof(SRepeatFilterOptions));
    if (*repeat_options == NULL)
        return BLASTERR_MEMORY;

    (*repeat_options)->database = strdup(kDefaultRepeatFilterDb);

    return 0;
}

Int2 SRepeatFilterOptionsResetDB(SRepeatFilterOptions* *repeat_options, const char* db)
{
    Int2 status=0;

    if (*repeat_options == NULL)
      status = SRepeatFilterOptionsNew(repeat_options);

    if (status)
      return status;

    sfree((*repeat_options)->database);
    (*repeat_options)->database = strdup(db);

    return status;
}

Int2 SWindowMaskerOptionsResetDB(SWindowMaskerOptions ** winmask_options, const char* db)
{
    Int2 status=0;
    
    if (*winmask_options == NULL)
        status = SWindowMaskerOptionsNew(winmask_options);
    
    if (status)
        return status;
    
    sfree((*winmask_options)->database);
    
    if (db) {
        (*winmask_options)->database = strdup(db);
    }
    
    return status;
}

SBlastFilterOptions* SBlastFilterOptionsFree(SBlastFilterOptions* filter_options)
{
    if (filter_options)
    {
        filter_options->dustOptions = 
            SDustOptionsFree(filter_options->dustOptions);
        filter_options->segOptions = 
            SSegOptionsFree(filter_options->segOptions);
        filter_options->repeatFilterOptions = 
            SRepeatFilterOptionsFree(filter_options->repeatFilterOptions);
        filter_options->windowMaskerOptions =
            SWindowMaskerOptionsFree(filter_options->windowMaskerOptions);
        sfree(filter_options);
    }

    return NULL;
}

Int2 SBlastFilterOptionsNew(SBlastFilterOptions* *filter_options,  EFilterOptions type)
{
    Int2 status = 0;

    if (filter_options)
    {
        *filter_options = (SBlastFilterOptions*) calloc(1, sizeof(SBlastFilterOptions));
        (*filter_options)->mask_at_hash = FALSE;
        if (type == eSeg)
          SSegOptionsNew(&((*filter_options)->segOptions)); 
        if (type == eDust || type == eDustRepeats)
          SDustOptionsNew(&((*filter_options)->dustOptions)); 
        if (type == eRepeats || type == eDustRepeats)
          SRepeatFilterOptionsNew(&((*filter_options)->repeatFilterOptions)); 
    }
    else
        status = 1;

    return status;
}


/** Merges together two sets of dust options, choosing the most non-default one.
 * 
 * @param opt1 first set to be merged [in]
 * @param opt2 second set to be merged [in]
 * @return the merged options.
 */
static SDustOptions* s_MergeDustOptions(const SDustOptions* opt1, const SDustOptions* opt2)
{
     SDustOptions* retval = NULL;

     if (!opt1 && !opt2)
         return NULL;

     SDustOptionsNew(&retval);

     if (opt1 && !opt2)
     {
           retval->level = opt1->level;
           retval->window = opt1->window;
           retval->linker = opt1->linker;
     }
     else if (!opt1 && opt2)
     {
           retval->level = opt2->level;
           retval->window = opt2->window;
           retval->linker = opt2->linker;
     }
     else 
     {
          retval->level = (opt1->level != kDustLevel) ? opt1->level : opt2->level;
          retval->window = (opt1->window != kDustWindow) ? opt1->window : opt2->window;
          retval->linker = (opt1->linker != kDustLinker) ? opt1->linker : opt2->linker;
     }

     return retval;
}


/** Merges together two sets of SEG options, choosing the most non-default one.
 * 
 * @param opt1 first set to be merged [in]
 * @param opt2 second set to be merged [in]
 * @return the merged options.
 */
static SSegOptions* s_MergeSegOptions(const SSegOptions* opt1, const SSegOptions* opt2)
{
    SSegOptions* retval = NULL;

    if (!opt1 && !opt2)
        return NULL;

    SSegOptionsNew(&retval);

    if (opt1 && !opt2)
    {
         retval->window = opt1->window;
         retval->locut = opt1->locut;
         retval->hicut = opt1->hicut;
    }
    else if (!opt1 && opt2)
    {
         retval->window = opt2->window;
         retval->locut = opt2->locut;
         retval->hicut = opt2->hicut;
    }
    else
    {
         retval->window = (opt1->window != kSegWindow) ? opt1->window : opt2->window;
         retval->locut = (opt1->locut != kSegLocut) ? opt1->locut : opt2->locut;
         retval->hicut = (opt1->hicut != kSegHicut) ? opt1->hicut : opt2->hicut;
    }
    return retval;
}

/** Merges together two sets of repeat filter options, choosing the most non-default one.
 * 
 * @param opt1 first set to be merged [in]
 * @param opt2 second set to be merged [in]
 * @return the merged options.
 */
static SRepeatFilterOptions* s_MergeRepeatOptions(const SRepeatFilterOptions* opt1, const SRepeatFilterOptions* opt2)
{
      SRepeatFilterOptions* retval = NULL;
  
      if (!opt1 && !opt2)
         return NULL;

      SRepeatFilterOptionsNew(&retval);

      if (opt1 && !opt2)
      {
           SRepeatFilterOptionsResetDB(&retval, opt1->database);
      }
      else if (!opt1 && opt2)
      {
           SRepeatFilterOptionsResetDB(&retval, opt2->database);
      }
      else 
      {  /* TODO : handle different db's. */
           SRepeatFilterOptionsResetDB(&retval, opt2->database); 
      }
      return retval;
}

/** Merges together two sets of window masker options, choosing the most non-default one.
 * 
 * @param opt1 first set to be merged [in]
 * @param opt2 second set to be merged [in]
 * @return the merged options.
 */
static SWindowMaskerOptions*
s_MergeWindowMaskerOptions(const SWindowMaskerOptions* opt1,
                           const SWindowMaskerOptions* opt2)
{
    SWindowMaskerOptions* retval = NULL;
    const SWindowMaskerOptions* src = NULL;
    Boolean have1 = FALSE, have2 = FALSE;
    
    have1 = opt1 && (opt1->database || opt1->taxid);
    have2 = opt2 && (opt2->database || opt2->taxid);
    
    if (! (have1 || have2))
        return NULL;
    
    if (have1 && ! have2) {
        src = opt1;
    } else if (! have1 && have2) {
        src = opt2;
    } else {
        // We have data structures with some kind of content, so
        // prefer structure 2 as repeat filter options do.
        src = opt2;
    }
    
    ASSERT(src);
    ASSERT(src->database || src->taxid);
    
    SWindowMaskerOptionsNew(&retval);
    SWindowMaskerOptionsResetDB(& retval, src->database);
    retval->taxid = src->taxid;
    
    return retval;
}

Int2 SBlastFilterOptionsMerge(SBlastFilterOptions** combined, const SBlastFilterOptions* opt1,
       const SBlastFilterOptions* opt2)
{
     SBlastFilterOptions* retval = NULL;
     Int2 status = 0;

     *combined = NULL;

     if (opt1 == NULL && opt2 == NULL)
         return 0;

     status = SBlastFilterOptionsNew(&retval, eEmpty);
     if (status != 0)
         return status;

     *combined = retval;

     if ((opt1 && opt1->mask_at_hash) || (opt2 && opt2->mask_at_hash))
         retval->mask_at_hash = TRUE;

     retval->dustOptions = 
         s_MergeDustOptions(opt1 ? opt1->dustOptions : NULL, opt2 ? opt2->dustOptions : NULL);
     retval->segOptions = 
         s_MergeSegOptions(opt1 ? opt1->segOptions : NULL, opt2 ? opt2->segOptions : NULL);
     retval->repeatFilterOptions = 
         s_MergeRepeatOptions(opt1 ? opt1->repeatFilterOptions : NULL, opt2 ? opt2->repeatFilterOptions : NULL);
     retval->windowMaskerOptions = 
         s_MergeWindowMaskerOptions(opt1 ? opt1->windowMaskerOptions : NULL, opt2 ? opt2->windowMaskerOptions : NULL);

     return 0;
}

Boolean SBlastFilterOptionsNoFiltering(const SBlastFilterOptions* filter_options)
{
       if (filter_options == NULL)
          return TRUE;
      
       return filter_options->dustOptions == NULL &&
           filter_options->segOptions == NULL &&
           filter_options->repeatFilterOptions == NULL &&
           filter_options->windowMaskerOptions == NULL;
}

Boolean SBlastFilterOptionsMaskAtHash(const SBlastFilterOptions* filter_options)
{
       if (filter_options == NULL)
          return FALSE;
      
       return filter_options->mask_at_hash;
}

Int2 SBlastFilterOptionsValidate(EBlastProgramType program_number, const SBlastFilterOptions* filter_options, Blast_Message* *blast_message)
{
       Int2 status = 0;

       if (filter_options == NULL)
       {
           Blast_MessageWrite(blast_message, eBlastSevWarning, kBlastMessageNoContext, 
              "SBlastFilterOptionsValidate: NULL filter_options");
           return BLASTERR_INVALIDPARAM;
       }

       if (filter_options->repeatFilterOptions)
       {
           if (program_number != eBlastTypeBlastn)
           {
               if (blast_message)
                  Blast_MessageWrite(blast_message, eBlastSevError, kBlastMessageNoContext,
                   "SBlastFilterOptionsValidate: Repeat filtering only supported with blastn");
               return  BLASTERR_OPTION_PROGRAM_INVALID;
           }
           if (filter_options->repeatFilterOptions->database == NULL ||
               strlen(filter_options->repeatFilterOptions->database) == 0)
           {
               if (blast_message)
                  Blast_MessageWrite(blast_message, eBlastSevError, kBlastMessageNoContext,
                   "SBlastFilterOptionsValidate: No repeat database specified for repeat filtering");
               return BLASTERR_INVALIDPARAM;
           }
       }

       if (filter_options->dustOptions)
       {
           if (program_number != eBlastTypeBlastn)
           {
               if (blast_message)
                  Blast_MessageWrite(blast_message, eBlastSevError, kBlastMessageNoContext,
                   "SBlastFilterOptionsValidate: Dust filtering only supported with blastn");
               return BLASTERR_OPTION_PROGRAM_INVALID;
           }
       }
  
       if (filter_options->segOptions)
       {
           if (program_number == eBlastTypeBlastn)
           {
               if (blast_message)
                  Blast_MessageWrite(blast_message, eBlastSevError, kBlastMessageNoContext,
                   "SBlastFilterOptionsValidate: SEG filtering is not supported with blastn");
               return BLASTERR_OPTION_PROGRAM_INVALID;
           }
       }

       return status;
}


QuerySetUpOptions*
BlastQuerySetUpOptionsFree(QuerySetUpOptions* options)

{
   if (options)
   {
       sfree(options->filter_string);
       options->filtering_options = SBlastFilterOptionsFree(options->filtering_options);
       sfree(options);
   }
   return NULL;
}

Int2
BlastQuerySetUpOptionsNew(QuerySetUpOptions* *options)
{
   Int2 status = 0;

   if (options == NULL)
      return BLASTERR_INVALIDPARAM;

   *options = (QuerySetUpOptions*) calloc(1, sizeof(QuerySetUpOptions));
   
   if (*options == NULL)
      return BLASTERR_MEMORY;

   (*options)->genetic_code = BLAST_GENETIC_CODE;

   /** @todo the code below should be deprecated */
   status = SBlastFilterOptionsNew(&((*options)->filtering_options), eEmpty);
   
   return status;
}

Int2 BLAST_FillQuerySetUpOptions(QuerySetUpOptions* options,
        EBlastProgramType program, const char *filter_string, Uint1 strand_option)
{
   Int2 status = 0;

   if (options == NULL)
      return BLASTERR_INVALIDPARAM;
   
   if (strand_option && 
       (program == eBlastTypeBlastn || program == eBlastTypePhiBlastn || 
        program == eBlastTypeBlastx || program == eBlastTypeTblastx)) {
      options->strand_option = strand_option;
   }

   if (filter_string) {
       /* Free whatever filter string has been set before. */
       sfree(options->filter_string);
       /* Free whatever filtering options have been set. */
       options->filtering_options =  SBlastFilterOptionsFree(options->filtering_options);
       /* Parse the filter_string for options, do not save the string. */
       status = BlastFilteringOptionsFromString(program, filter_string, 
          &options->filtering_options, NULL);
   }
   return status;
}

BlastInitialWordOptions*
BlastInitialWordOptionsFree(BlastInitialWordOptions* options)

{

	sfree(options);

	return NULL;
}


Int2
BlastInitialWordOptionsNew(EBlastProgramType program, 
   BlastInitialWordOptions* *options)
{
   *options = 
      (BlastInitialWordOptions*) calloc(1, sizeof(BlastInitialWordOptions));
   if (*options == NULL)
      return BLASTERR_MEMORY;

   if (program != eBlastTypeBlastn &&
       program != eBlastTypePhiBlastn) {	/* protein-protein options. */
      (*options)->window_size = BLAST_WINDOW_SIZE_PROT;
      (*options)->x_dropoff = BLAST_UNGAPPED_X_DROPOFF_PROT;
      (*options)->gap_trigger = BLAST_GAP_TRIGGER_PROT;
   } else {
      (*options)->window_size = BLAST_WINDOW_SIZE_NUCL;
      (*options)->scan_range =  BLAST_SCAN_RANGE_NUCL;
      (*options)->gap_trigger = BLAST_GAP_TRIGGER_NUCL;
      (*options)->x_dropoff = BLAST_UNGAPPED_X_DROPOFF_NUCL;
   }

   (*options)->program_number = program;

   return 0;
}


Int2
BlastInitialWordOptionsValidate(EBlastProgramType program_number,
   const BlastInitialWordOptions* options, 
   Blast_Message* *blast_msg)
{

   ASSERT(options);

   /* PHI-BLAST has no ungapped extension phase.  Megablast may not have it,
    but generally does now. */
   if (program_number != eBlastTypeBlastn  &&
       (!Blast_ProgramIsPhiBlast(program_number)) &&
       options->x_dropoff <= 0.0)
   {
      Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                            "x_dropoff must be greater than zero");
         return BLASTERR_OPTION_VALUE_INVALID;
   }

   if (program_number == eBlastTypeBlastn && 
       options->scan_range && !options->window_size)
   {
      Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                            "off_diagonal_range is only useful in 2-hit algorithm");
         return BLASTERR_OPTION_VALUE_INVALID;
   }
                
   
   return 0;
}


Int2
BLAST_FillInitialWordOptions(BlastInitialWordOptions* options, 
                EBlastProgramType program, Int4 window_size, 
                double xdrop_ungapped)
{
   if (!options)
      return BLASTERR_INVALIDPARAM;

   if (window_size != 0)
      options->window_size = window_size;
   if (xdrop_ungapped != 0)
      options->x_dropoff = xdrop_ungapped;

   return 0;
}

BlastExtensionOptions*
BlastExtensionOptionsFree(BlastExtensionOptions* options)

{

	sfree(options);

	return NULL;
}

BlastScoringOptions*
BlastScoringOptionsFree(BlastScoringOptions* options)

{
	if (options == NULL)
		return NULL;

	sfree(options->matrix);
   sfree(options->matrix_path);
	sfree(options);

	return NULL;
}

Int2 
BlastScoringOptionsNew(EBlastProgramType program_number, BlastScoringOptions* *options)
{
   *options = (BlastScoringOptions*) calloc(1, sizeof(BlastScoringOptions));

   if (*options == NULL)
      return BLASTERR_INVALIDPARAM;
   
   if (program_number != eBlastTypeBlastn &&
       program_number != eBlastTypePhiBlastn) {	/* protein-protein options. */
      (*options)->shift_pen = INT2_MAX;
      (*options)->is_ooframe = FALSE;
      (*options)->gap_open = BLAST_GAP_OPEN_PROT;
      (*options)->gap_extend = BLAST_GAP_EXTN_PROT;
      (*options)->matrix = strdup(BLAST_DEFAULT_MATRIX);
   } else {	/* nucleotide-nucleotide options. */
      (*options)->penalty = BLAST_PENALTY;
      (*options)->reward = BLAST_REWARD;
      /* This is correct except when greedy extension is used. In that case 
         these values would have to be reset. */
      (*options)->gap_open = BLAST_GAP_OPEN_NUCL;
      (*options)->gap_extend = BLAST_GAP_EXTN_NUCL;
   }
   if (program_number != eBlastTypeTblastx) {
       (*options)->gapped_calculation = TRUE;
   }
   (*options)->program_number = program_number;
   /* By default cross_match-like complexity adjusted scoring is 
      turned off.  RMBlastN is currently the only program to use this. -RMH */  
   (*options)->complexity_adjusted_scoring = FALSE;
   
   return 0;
}

Int2 
BLAST_FillScoringOptions(BlastScoringOptions* options, 
   EBlastProgramType program_number, Boolean greedy_extension, Int4 penalty, Int4 reward, 
   const char *matrix, Int4 gap_open, Int4 gap_extend)
{
   if (!options)
      return BLASTERR_INVALIDPARAM;

   if (program_number != eBlastTypeBlastn &&
       program_number != eBlastTypePhiBlastn) {	/* protein-protein options. */
      /* If matrix name is not provided, keep the default "BLOSUM62" value filled in 
         BlastScoringOptionsNew, otherwise reset it. */
      if (matrix)
          BlastScoringOptionsSetMatrix(options, matrix);
   } else {	/* nucleotide-nucleotide options. */
      if (penalty)
         options->penalty = penalty;
      if (reward)
         options->reward = reward;

      if (greedy_extension) {
         options->gap_open = BLAST_GAP_OPEN_MEGABLAST;
         options->gap_extend = BLAST_GAP_EXTN_MEGABLAST;
      }	else {
         options->gap_open = BLAST_GAP_OPEN_NUCL;
         options->gap_extend = BLAST_GAP_EXTN_NUCL;
      }
   }
   if (gap_open >= 0)
      options->gap_open = gap_open;
   if (gap_extend >= 0)
      options->gap_extend = gap_extend;

   options->program_number = program_number;

   return 0;
}

Int2 
BlastScoringOptionsValidate(EBlastProgramType program_number, 
   const BlastScoringOptions* options, Blast_Message* *blast_msg)

{
	if (options == NULL)
		return BLASTERR_INVALIDPARAM;

        if (program_number == eBlastTypeTblastx && options->gapped_calculation)
        {
            Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
               "Gapped search is not allowed for tblastx");
		return BLASTERR_OPTION_PROGRAM_INVALID;
        }

	if (program_number == eBlastTypeBlastn || program_number == eBlastTypePhiBlastn)
	{
           // A penalty/reward of 0/0 is a signal that this is rmblastn
           // which allows specification of penalties as positive integers.
           if ( ! ( options->penalty == 0 && options->reward == 0 ) )
           {
		if (options->penalty >= 0)
		{
			Blast_MessageWrite(blast_msg, eBlastSevWarning, kBlastMessageNoContext,
                            "BLASTN penalty must be negative");
			return BLASTERR_OPTION_VALUE_INVALID;
		}

                if (options->gapped_calculation && !BLAST_CheckRewardPenaltyScores(options->reward, options->penalty))
                {
			Blast_MessageWrite(blast_msg, eBlastSevWarning, kBlastMessageNoContext,
                            "BLASTN reward/penalty combination not supported for gapped search");
			return BLASTERR_OPTION_VALUE_INVALID;
                }
             }

             if (options->gapped_calculation && options->gap_open > 0 && options->gap_extend == 0) 
             {
                     Blast_MessageWrite(blast_msg, eBlastSevWarning, kBlastMessageNoContext,
                        "BLASTN gap extension penalty cannot be 0");
                     return BLASTERR_OPTION_VALUE_INVALID;
             }
	}
	else
	{
                if (options->gapped_calculation && !Blast_ProgramIsRpsBlast(program_number))
                {
                    Int2 status=0;
                    if ((status=Blast_KarlinBlkGappedLoadFromTables(NULL, options->gap_open,
                     options->gap_extend, options->matrix)) != 0)
                     {
			if (status == 1)
			{
				char* buffer;

				buffer = BLAST_PrintMatrixMessage(options->matrix); 
                                Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext, buffer);
				sfree(buffer);
				return BLASTERR_OPTION_VALUE_INVALID;
				
			}
			else if (status == 2)
			{
				char* buffer;

				buffer = BLAST_PrintAllowedValues(options->matrix, 
                        options->gap_open, options->gap_extend);
                                Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext, buffer);
				sfree(buffer);
				return BLASTERR_OPTION_VALUE_INVALID;
			}
                    }
	       }
	}

	if (program_number != eBlastTypeBlastx && program_number != eBlastTypeTblastn && options->is_ooframe)
	{
            Blast_MessageWrite(blast_msg, eBlastSevWarning, kBlastMessageNoContext,
               "Out-of-frame only permitted for blastx and tblastn");
            return  BLASTERR_OPTION_PROGRAM_INVALID;
	}

	return 0;
}

Int2 
BlastScoringOptionsDup(BlastScoringOptions* *new_opt, const BlastScoringOptions* old_opt)
{
    if (old_opt == NULL || new_opt == NULL)
       return BLASTERR_INVALIDPARAM;

    *new_opt = (BlastScoringOptions*) BlastMemDup(old_opt, sizeof(BlastScoringOptions));
    if (*new_opt == NULL)
       return BLASTERR_MEMORY;

    if (old_opt->matrix)
       (*new_opt)->matrix = strdup(old_opt->matrix);

    if (old_opt->matrix_path)
       (*new_opt)->matrix_path = strdup(old_opt->matrix_path);

    return 0;
}

Int2 BlastScoringOptionsSetMatrix(BlastScoringOptions* opts,
                                  const char* matrix_name)
{
    Uint4 i;

    if (matrix_name) {
        sfree(opts->matrix);
        opts->matrix = strdup(matrix_name);
        /* Make it all upper case */
        for (i=0; i<strlen(opts->matrix); ++i)
            opts->matrix[i] = toupper((unsigned char) opts->matrix[i]);
    }
    return 0;
}

BlastEffectiveLengthsOptions*
BlastEffectiveLengthsOptionsFree(BlastEffectiveLengthsOptions* options)

{
   if (options == NULL)
      return NULL;

   sfree(options->searchsp_eff);
   sfree(options);
   return NULL;
}


Int2 
BlastEffectiveLengthsOptionsNew(BlastEffectiveLengthsOptions* *options)

{
    if (options == NULL) {
        return BLASTERR_INVALIDPARAM;
    }

    *options = (BlastEffectiveLengthsOptions*)
       calloc(1, sizeof(BlastEffectiveLengthsOptions));
 
    if (*options == NULL)
       return BLASTERR_MEMORY;
    
    return 0;
}

Boolean
BlastEffectiveLengthsOptions_IsSearchSpaceSet(const
                                              BlastEffectiveLengthsOptions*
                                              options)
{
    int i;
    if ( !options || options->searchsp_eff == NULL) {
        return FALSE;
    }

    for (i = 0; i < options->num_searchspaces; i++) {
        if (options->searchsp_eff[i] != 0) {
            return TRUE;
        }
    }
    return FALSE;
}

Int2 
BLAST_FillEffectiveLengthsOptions(BlastEffectiveLengthsOptions* options, 
   Int4 dbseq_num, Int8 db_length, Int8* searchsp_eff, Int4 num_searchsp)
{
   Int4 index;
   if (!options)
      return BLASTERR_INVALIDPARAM;

   if (num_searchsp > options->num_searchspaces) {
       options->num_searchspaces = num_searchsp;
       options->searchsp_eff = (Int8 *)realloc(options->searchsp_eff,
                                               num_searchsp * sizeof(Int8));
       if (options->searchsp_eff == NULL)
           return BLASTERR_MEMORY;
   }

   for (index = 0; index < options->num_searchspaces; index++)
      options->searchsp_eff[index] = searchsp_eff[index];

   options->dbseq_num = dbseq_num;
   options->db_length = db_length;

   return 0;
}

LookupTableOptions*
LookupTableOptionsFree(LookupTableOptions* options)

{

      if (options == NULL)
          return NULL;

      sfree(options->phi_pattern);
   
	sfree(options);
	return NULL;
}

Int2 
LookupTableOptionsNew(EBlastProgramType program_number, LookupTableOptions* *options)
{
   *options = (LookupTableOptions*) calloc(1, sizeof(LookupTableOptions));
   
   if (*options == NULL)
      return BLASTERR_INVALIDPARAM;
   
   switch (program_number) {
   case eBlastTypeBlastn:
       /* Blastn default is megablast. */
       (*options)->word_size = BLAST_WORDSIZE_MEGABLAST;
       (*options)->lut_type = eMBLookupTable;
       break;
   case eBlastTypeRpsBlast: case eBlastTypeRpsTblastn:
       (*options)->word_size = BLAST_WORDSIZE_PROT;
       (*options)->lut_type = eRPSLookupTable;
       
       if (program_number == eBlastTypeRpsBlast)
           (*options)->threshold = BLAST_WORD_THRESHOLD_BLASTP;
       else 
           (*options)->threshold = BLAST_WORD_THRESHOLD_TBLASTN;
       break;
   case eBlastTypePhiBlastn:
       (*options)->lut_type = ePhiNaLookupTable;
       break;
   case eBlastTypePhiBlastp:
       (*options)->lut_type = ePhiLookupTable;
       break;
   default:
       (*options)->word_size = BLAST_WORDSIZE_PROT;
       (*options)->lut_type = eAaLookupTable;
       
       if (program_number == eBlastTypeBlastp)
           (*options)->threshold = BLAST_WORD_THRESHOLD_BLASTP;
       else if (program_number == eBlastTypeBlastx)
           (*options)->threshold = BLAST_WORD_THRESHOLD_BLASTX;
       else if (program_number == eBlastTypeTblastn)
           (*options)->threshold = BLAST_WORD_THRESHOLD_TBLASTN;
       else if (program_number == eBlastTypeTblastx)
           (*options)->threshold = BLAST_WORD_THRESHOLD_TBLASTX;
       break;
   }

   (*options)->program_number = program_number;

   return 0;
}

Int2 
BLAST_FillLookupTableOptions(LookupTableOptions* options, 
   EBlastProgramType program_number, Boolean is_megablast, 
   double threshold, Int4 word_size)
{
   if (!options)
      return BLASTERR_INVALIDPARAM;

   if (program_number == eBlastTypeBlastn) {
      if (is_megablast)	{
         options->lut_type = eMBLookupTable;
         options->word_size = BLAST_WORDSIZE_MEGABLAST;
      }	else {
         options->lut_type = eNaLookupTable;
         options->word_size = BLAST_WORDSIZE_NUCL;
      }
   } else {
      options->lut_type = eAaLookupTable;
   }

   /* if the supplied threshold is negative, disable neighboring words */
   if (threshold < 0)
      options->threshold = 0;

   /* if the supplied threshold is > 0, use it otherwise, use the default */
   if (threshold > 0)
      options->threshold = threshold;

   if (Blast_ProgramIsRpsBlast(program_number))
      options->lut_type = eRPSLookupTable;
   if (word_size)
      options->word_size = word_size;
   if ((program_number == eBlastTypeTblastn ||
        program_number == eBlastTypeBlastp ||
        program_number == eBlastTypeBlastx) && 
       word_size > 5)
       options->lut_type = eCompressedAaLookupTable;

   return 0;
}

Int2 BLAST_GetSuggestedThreshold(EBlastProgramType program_number, const char* matrixName, double* threshold)
{

    const double kB62_threshold = 11;

    if (program_number == eBlastTypeBlastn)
      return 0;

    if (matrixName == NULL)
      return BLASTERR_INVALIDPARAM;

    if(strcasecmp(matrixName, "BLOSUM62") == 0)
        *threshold = kB62_threshold;
    else if(strcasecmp(matrixName, "BLOSUM45") == 0)
        *threshold = 14;
    else if(strcasecmp(matrixName, "BLOSUM62_20") == 0)
        *threshold = 100;
    else if(strcasecmp(matrixName, "BLOSUM80") == 0)
        *threshold = 12;
    else if(strcasecmp(matrixName, "PAM30") == 0)
        *threshold = 16;
    else if(strcasecmp(matrixName, "PAM70") == 0)
        *threshold = 14;
    else
        *threshold = kB62_threshold;

    if (Blast_SubjectIsTranslated(program_number) == TRUE)
        *threshold += 2;  /* Covers tblastn, tblastx, psi-tblastn rpstblastn. */
    else if (Blast_QueryIsTranslated(program_number) == TRUE)
        *threshold += 1;

    return 0;
}

Int2 BLAST_GetSuggestedWindowSize(EBlastProgramType program_number, const char* matrixName, Int4* window_size)
{
    const Int4 kB62_windowsize = 40;

    if (program_number == eBlastTypeBlastn)
      return 0;

    if (matrixName == NULL)
      return BLASTERR_INVALIDPARAM;

    if(strcasecmp(matrixName, "BLOSUM62") == 0)
        *window_size = kB62_windowsize;
    else if(strcasecmp(matrixName, "BLOSUM45") == 0)
        *window_size = 60;
    else if(strcasecmp(matrixName, "BLOSUM80") == 0)
        *window_size = 25;
    else if(strcasecmp(matrixName, "PAM30") == 0)
        *window_size = 15;
    else if(strcasecmp(matrixName, "PAM70") == 0)
        *window_size = 20;
    else
        *window_size = kB62_windowsize;

    return 0;
}

/** Validate options for the discontiguous word megablast
 * Word size must be 11 or 12; template length 16, 18 or 21; 
 * template type 0, 1 or 2.
 * @param word_size Word size option [in]
 * @param template_length Discontiguous template length [in]
 * @param template_type Discontiguous template type [in]
 * @param blast_msg Used for storing error messages [in][out]
 * @return TRUE if options combination valid.
 */
static Boolean 
s_DiscWordOptionsValidate(Int4 word_size, Uint1 template_length,
                          Uint1 template_type,
                          Blast_Message** blast_msg)
{
   if (template_length == 0)
      return TRUE;


   if (word_size != 11 && word_size != 12) {
      Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                         "Invalid discontiguous template parameters: word "
                         "size must be either 11 or 12");
      return FALSE;
   }

   if (template_length != 16 && template_length != 18 && 
       template_length != 21) {
      Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                         "Invalid discontiguous template parameters: "
                         "template length must be 16, 18, or 21");
      return FALSE;
   }

   if (template_type > 2) {
     /* should never fail coming from the C++ APIs as we represent these as
      * strings */
      Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                         "Invalid discontiguous template parameters: "
                         "template type must be 0, 1, or 2");
      return FALSE;
   }

   return TRUE;
}

Int2 
LookupTableOptionsValidate(EBlastProgramType program_number, 
   const LookupTableOptions* options, Blast_Message* *blast_msg)

{
   const Boolean kPhiBlast = Blast_ProgramIsPhiBlast(program_number);

    if (options == NULL)
        return BLASTERR_INVALIDPARAM;

    if (options->phi_pattern && !kPhiBlast) {
        Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
            "PHI pattern can be specified only for blastp and blastn");
        return BLASTERR_OPTION_PROGRAM_INVALID;
    }

    /* For PHI BLAST, the subsequent word size tests are not needed. */
    if (kPhiBlast)
        return 0;

    if (program_number != eBlastTypeBlastn && 
        (!Blast_ProgramIsRpsBlast(program_number)) &&
        options->threshold <= 0)
    {
        Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                         "Non-zero threshold required");
        return BLASTERR_OPTION_VALUE_INVALID;
    }

    if (options->word_size <= 0)
    {
        if ( !Blast_ProgramIsRpsBlast(program_number)) {
            Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                                     "Word-size must be greater than zero");
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    } else if (program_number == eBlastTypeBlastn && options->word_size < 4)
    {
        Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext, 
                  "Word-size must be 4 or greater for nucleotide comparison");
        return BLASTERR_OPTION_VALUE_INVALID;
    } else if (program_number != eBlastTypeBlastn && options->word_size > 5)
    {
        if (program_number == eBlastTypeBlastp ||
            program_number == eBlastTypeTblastn ||
            program_number == eBlastTypeBlastx)
        {
            if (options->word_size > 7) {
                Blast_MessageWrite(blast_msg, eBlastSevError, 
                                   kBlastMessageNoContext,
                                   "Word-size must be less than "
                                   "8 for a tblastn, blastp or blastx search");
                return BLASTERR_OPTION_VALUE_INVALID;
            }
        }
        else {
            Blast_MessageWrite(blast_msg, eBlastSevError, 
                               kBlastMessageNoContext,
                               "Word-size must be less "
                               "than 6 for protein comparison");
            return BLASTERR_OPTION_VALUE_INVALID;
        }
    }

    if (program_number != eBlastTypeBlastn && 
       options->lut_type == eMBLookupTable)
    {
        Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                         "Megablast lookup table only supported with blastn");
        return BLASTERR_OPTION_PROGRAM_INVALID;
    }

    if (program_number == eBlastTypeBlastp ||
        program_number == eBlastTypeTblastn ||
        program_number == eBlastTypeBlastx)
    {
        if (options->word_size > 5 &&
            options->lut_type != eCompressedAaLookupTable) {
           Blast_MessageWrite(blast_msg, eBlastSevError,
                              kBlastMessageNoContext,
                              "Blastp, Blastx or Tblastn with word size"
                              " > 5 requires a "
                              "compressed alphabet lookup table");
           return BLASTERR_OPTION_VALUE_INVALID;
        }
        else if (options->lut_type == eCompressedAaLookupTable &&
                 options->word_size != 6 && options->word_size != 7) {
           Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                         "Compressed alphabet lookup table requires "
                         "word size 6 or 7");
           return BLASTERR_OPTION_VALUE_INVALID;
        }
    }

   if (program_number == eBlastTypeBlastn && options->mb_template_length > 0) {
      if (!s_DiscWordOptionsValidate(options->word_size,
              options->mb_template_length, 
              options->mb_template_type,
              blast_msg)) {
         return BLASTERR_OPTION_VALUE_INVALID;
      } else if (options->lut_type != eMBLookupTable) {
         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
            "Invalid lookup table type for discontiguous Mega BLAST");
         return BLASTERR_OPTION_VALUE_INVALID;
      } 
   }
    return 0;
}

BlastHitSavingOptions*
BlastHitSavingOptionsFree(BlastHitSavingOptions* options)

{
    if (options) {
        options->hsp_filt_opt = BlastHSPFilteringOptionsFree(options->hsp_filt_opt);
    }
    sfree(options);
    return NULL;
}


Int2 BlastHitSavingOptionsNew(EBlastProgramType program_number, 
        BlastHitSavingOptions** options,
        Boolean gapped_calculation)
{
   *options = (BlastHitSavingOptions*) calloc(1, sizeof(BlastHitSavingOptions));
   
   if (*options == NULL)
      return BLASTERR_INVALIDPARAM;

   (*options)->hitlist_size = BLAST_HITLIST_SIZE;
   (*options)->expect_value = BLAST_EXPECT_VALUE;
   (*options)->program_number = program_number;

   // Initialize mask_level parameter -RMH-
   (*options)->mask_level = 101;

   /* By default, sum statistics is used for all translated searches 
    * (except RPS BLAST), and for all ungapped searches.
    */
   if (program_number == eBlastTypeRpsTblastn) {
	   (*options)->do_sum_stats = FALSE;
   } else if (!gapped_calculation ||
	   Blast_QueryIsTranslated(program_number) ||
	   Blast_SubjectIsTranslated(program_number)) {
       (*options)->do_sum_stats = TRUE;
   } else {
       (*options)->do_sum_stats = FALSE;
   }

   (*options)->hsp_filt_opt = NULL;

   return 0;

}

Int2
BLAST_FillHitSavingOptions(BlastHitSavingOptions* options, 
                           double evalue, Int4 hitlist_size,
                           Boolean is_gapped, Int4 culling_limit,
                           Int4 min_diag_separation)
{
   if (!options)
      return BLASTERR_INVALIDPARAM;

   if (hitlist_size)
      options->hitlist_size = hitlist_size;
   if (evalue)
      options->expect_value = evalue;
   if (min_diag_separation)
      options->min_diag_separation = min_diag_separation;
   options->culling_limit = culling_limit;
   options->hsp_filt_opt = NULL;

   return 0;

}

Int2 PSIBlastOptionsNew(PSIBlastOptions** psi_options)
{
   PSIBlastOptions* options = NULL;

   if ( !psi_options )
      return BLASTERR_INVALIDPARAM;

   options = (PSIBlastOptions*)calloc(1, sizeof(PSIBlastOptions));
   if ( !options ) 
       return BLASTERR_MEMORY;

   *psi_options = options;
   options->inclusion_ethresh = PSI_INCLUSION_ETHRESH;
   options->pseudo_count = PSI_PSEUDO_COUNT_CONST;
   options->use_best_alignment = TRUE;

   options->nsg_compatibility_mode = FALSE;
   options->impala_scaling_factor = kPSSM_NoImpalaScaling;
   options->ignore_unaligned_positions = FALSE;
   
   return 0;
}

Int2 PSIBlastOptionsValidate(const PSIBlastOptions* psi_options,
                             Blast_Message** blast_msg)
{
    Int2 retval = 1;    /* assume failure */

    if ( !psi_options ) {
        return retval;
    }

    if (psi_options->pseudo_count < 0) {
        Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                           "Pseudo count must be greater than or equal to 0");
        return retval;
    }

    if (psi_options->inclusion_ethresh <= 0.0) {
        Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
                           "Inclusion threshold must be greater than 0");
        return retval;
    }

    retval = 0;
    return retval;
}

PSIBlastOptions* PSIBlastOptionsFree(PSIBlastOptions* psi_options)
{
   sfree(psi_options);
   return NULL;
}

Int2 BlastDatabaseOptionsNew(BlastDatabaseOptions** db_options)
{
   BlastDatabaseOptions* options = NULL;

   if ( !db_options ) {
       return BLASTERR_INVALIDPARAM;
   }

   options = (BlastDatabaseOptions*) calloc(1, sizeof(BlastDatabaseOptions));
   if ( !options ) {
       return  BLASTERR_MEMORY;
   }

   options->genetic_code = BLAST_GENETIC_CODE;
   *db_options = options;

   return 0;
}

BlastDatabaseOptions* 
BlastDatabaseOptionsFree(BlastDatabaseOptions* db_options)
{

   if (db_options == NULL)
      return NULL;

   sfree(db_options);
   return NULL;
}

BlastHSPBestHitOptions* BlastHSPBestHitOptionsNew(double overhang, double score_edge)
{
    BlastHSPBestHitOptions* retval = 
        (BlastHSPBestHitOptions*) calloc(1, sizeof(BlastHSPBestHitOptions));
    retval->overhang = overhang;
    retval->score_edge = score_edge;
    return retval;
}

BlastHSPBestHitOptions* BlastHSPBestHitOptionsFree(BlastHSPBestHitOptions* opt)
{
    if ( !opt ) {
        return NULL;
    }
    sfree(opt);
    return NULL;
}

BlastHSPCullingOptions* BlastHSPCullingOptionsNew(int max)
{
    BlastHSPCullingOptions* retval = 
        (BlastHSPCullingOptions*) calloc(1, sizeof(BlastHSPCullingOptions));
    retval->max_hits = max;
    return retval;
}

Int2
BlastHSPCullingOptionsValidate(const BlastHSPFilteringOptions* opts)
{
    Int2 retval = 0;
    BlastHSPCullingOptions* culling_opts = opts->culling_opts;
    if (!culling_opts)
       return retval;

    if (culling_opts->max_hits < 0)
       return -1;

    return retval;
}

BlastHSPCullingOptions* 
BlastHSPCullingOptionsFree(BlastHSPCullingOptions* culling_opts)
{
   if (!culling_opts)
    return NULL;

   sfree(culling_opts);
   return NULL;
}


BlastHSPFilteringOptions* BlastHSPFilteringOptionsNew()
{
    return (BlastHSPFilteringOptions*)calloc(1,
                                             sizeof(BlastHSPFilteringOptions));
}

Int2
BlastHSPFilteringOptions_AddBestHit(BlastHSPFilteringOptions* filt_opts,
                                    BlastHSPBestHitOptions** best_hit,
                                    EBlastStage stage)
{
    if ( filt_opts == NULL || best_hit == NULL || *best_hit == NULL) {
        return 1;
    }

    filt_opts->best_hit = *best_hit;
    *best_hit = NULL;
    filt_opts->best_hit_stage = stage;

    return 0;
}

Int2
BlastHSPFilteringOptions_AddCulling(BlastHSPFilteringOptions* filt_opts,
                                    BlastHSPCullingOptions** culling,
                                    EBlastStage stage)
{
    if ( filt_opts == NULL || culling == NULL || *culling == NULL) {
        return 1;
    }

    filt_opts->culling_opts = *culling;
    *culling = NULL;
    filt_opts->culling_stage = stage;

    return 0;
}

BlastHSPFilteringOptions*
BlastHSPFilteringOptionsFree(BlastHSPFilteringOptions* opts)
{
    if ( !opts ) {
        return NULL;
    }
    opts->best_hit = BlastHSPBestHitOptionsFree(opts->best_hit);
    opts->culling_opts = BlastHSPCullingOptionsFree(opts->culling_opts);
    sfree(opts);
    return opts;
}

