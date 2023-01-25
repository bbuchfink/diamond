/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: WaveFront-Alignment module for debugging and collect stats
 */

#include "../utils/commons.h"
#include "wavefront_debug.h"
#include "wavefront_align.h"
#include "wavefront_compute.h"

/*
 * Checks
 */
bool wavefront_check_alignment(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  const char* const pattern = wf_aligner->pattern;
  const int pattern_length = wf_aligner->pattern_length;
  const char* const text = wf_aligner->text;
  const int text_length = wf_aligner->text_length;
  // Custom function to compare sequences
  alignment_match_funct_t match_funct = wf_aligner->match_funct;
  void* match_funct_arguments = wf_aligner->match_funct_arguments;
  // CIGAR
  cigar_t* const cigar = wf_aligner->cigar;
  char* const operations = cigar->operations;
  const int begin_offset = cigar->begin_offset;
  const int end_offset = cigar->end_offset;
  // Traverse CIGAR
  bool alignment_correct = true;
  int pattern_pos=0, text_pos=0, i;
  for (i=begin_offset;i<end_offset;++i) {
    switch (operations[i]) {
      case 'M': {
        // Check match
        const bool is_match = (match_funct!=NULL) ?
            match_funct(pattern_pos,text_pos,match_funct_arguments) :
            pattern[pattern_pos] == text[text_pos];
        if (!is_match) {
          fprintf(stream,"[WFA::Check] Alignment not matching (pattern[%d]=%c != text[%d]=%c)\n",
              pattern_pos,pattern[pattern_pos],text_pos,text[text_pos]);
          alignment_correct = false;
          break;
        }
        ++pattern_pos;
        ++text_pos;
        break;
      }
      case 'X': {
        // Check mismatch
        const bool is_match = (match_funct!=NULL) ?
            match_funct(pattern_pos,text_pos,match_funct_arguments) :
            pattern[pattern_pos] == text[text_pos];
        if (is_match) {
          fprintf(stream,"[WFA::Check] Alignment not mismatching (pattern[%d]=%c == text[%d]=%c)\n",
              pattern_pos,pattern[pattern_pos],text_pos,text[text_pos]);
          alignment_correct = false;
          break;
        }
        ++pattern_pos;
        ++text_pos;
        break;
      }
      case 'I':
        ++text_pos;
        break;
      case 'D':
        ++pattern_pos;
        break;
      default:
        fprintf(stream,"[WFA::Check] Unknown edit operation '%c'\n",operations[i]);
        exit(1);
        break;
    }
  }
  // Check alignment length
  if (pattern_pos != pattern_length) {
    fprintf(stream,
        "[WFA::Check] Alignment incorrect length (pattern-aligned=%d,pattern-length=%d)\n",
        pattern_pos,pattern_length);
    alignment_correct = false;
  }
  if (text_pos != text_length) {
    fprintf(stream,
        "[WFA::Check] Alignment incorrect length (text-aligned=%d,text-length=%d)\n",
        text_pos,text_length);
    alignment_correct = false;
  }
  // Return
  return alignment_correct;
}
/*
 * Reporting
 */
void wavefront_report_lite(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  const char* const pattern = wf_aligner->pattern;
  const int pattern_length = wf_aligner->pattern_length;
  const char* const text = wf_aligner->text;
  const int text_length = wf_aligner->text_length;
  const int status = wf_aligner->align_status.status;
  const uint64_t memory_used = wf_aligner->align_status.memory_used;
  // Banner
  fprintf(stream,"[WFA::Debug]");
  // Sequences
  const int score = wavefront_compute_classic_score(
      wf_aligner,wf_aligner->pattern_length,
      wf_aligner->text_length,wf_aligner->cigar->score);
  fprintf(stream,"\t%d",score);
  fprintf(stream,"\t%d\t%d",pattern_length,text_length);
  fprintf(stream,"\t%s",(status==0) ? "OK" : "FAIL");
  fprintf(stream,"\t%2.3f",TIMER_GET_TOTAL_MS(&wf_aligner->system.timer));
  fprintf(stream,"\t%luMB\t",CONVERT_B_TO_MB(memory_used));
  fprintf(stream,"[");
  wavefront_aligner_print_type(stream,wf_aligner);
  fprintf(stream,",");
  wavefront_aligner_print_scope(stream,wf_aligner);
  fprintf(stream,",");
  wavefront_penalties_print(stream,&wf_aligner->penalties);
  fprintf(stream,"]\t");
  cigar_print(stream,wf_aligner->cigar,true);
  if (wf_aligner->match_funct != NULL) {
    fprintf(stream,"\t-\t-");
  } else {
    fprintf(stream,"\t%.*s\t%.*s",pattern_length,pattern,text_length,text);
  }
  fprintf(stream,"\n");
}
void wavefront_report_verbose_begin(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  // Input sequences
  fprintf(stream,"[WFA::Report::Begin] [");
  wavefront_aligner_print_type(stream,wf_aligner);
  fprintf(stream,"]-Alignment (obj=%p)\n",wf_aligner);
  if (wf_aligner->match_funct != NULL) {
    fprintf(stream,"[WFA::Report]\tPattern\t%d\tcustom-funct()\n",pattern_length);
    fprintf(stream,"[WFA::Report]\tText\t%d\tcustom-funct()\n",text_length);
  } else {
    fprintf(stream,"[WFA::Report]\tPattern\t%d\t%.*s\n",pattern_length,pattern_length,pattern);
    fprintf(stream,"[WFA::Report]\tText\t%d\t%.*s\n",text_length,text_length,text);
  }
  // Alignment scope/form
  fprintf(stream,"[WFA::Report]\tScope=");
  wavefront_aligner_print_scope(stream,wf_aligner);
  fprintf(stream," Max-score=%d",
      wf_aligner->system.max_alignment_score);
  // Penalties
  fprintf(stream," Penalties=");
  wavefront_penalties_print(stream,&wf_aligner->penalties);
  // Heuristic
  fprintf(stream," Heuristic=");
  wavefront_heuristic_print(stream,&wf_aligner->heuristic);
  // Memory mode
  fprintf(stream," Memory.mode=(%d,%luMB,%luMB,%luMB)\n",
      wf_aligner->memory_mode,
      CONVERT_B_TO_MB(wf_aligner->system.max_memory_compact),
      CONVERT_B_TO_MB(wf_aligner->system.max_memory_resident),
      CONVERT_B_TO_MB(wf_aligner->system.max_memory_abort));
}
void wavefront_report_verbose_end(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner) {
  // Finish report
  fprintf(stream,"[WFA::Report::End]\tFinish.status=%d",wf_aligner->align_status.status);
  fprintf(stream," Time.taken=");
  timer_print_total(stream,&wf_aligner->system.timer);
  fprintf(stream," Memory.used=%luMB",
      CONVERT_B_TO_MB(wf_aligner->align_status.memory_used));
  fprintf(stream," WFA.components=(wfs=%d,maxlo=%d,maxhi=%d)",
      wf_aligner->wf_components.num_wavefronts,
      wf_aligner->wf_components.historic_min_lo,
      wf_aligner->wf_components.historic_max_hi);
  const int score = wavefront_compute_classic_score(
      wf_aligner,wf_aligner->pattern_length,
      wf_aligner->text_length,wf_aligner->cigar->score);
  fprintf(stream," WFA.score=%d",score);
  fprintf(stream," WFA.cigar=");
  cigar_print(stream,wf_aligner->cigar,true);
  fprintf(stream,"\n");
}
/*
 * Debug
 */
void wavefront_debug_prologue(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  // Check verbose level
  if (wf_aligner->system.verbose >= 1) {
    timer_start(&wf_aligner->system.timer);
    if (wf_aligner->system.verbose >= 4) {
      wavefront_report_verbose_begin(stderr,wf_aligner,
          pattern,pattern_length,text,text_length);
    }
  }
}
void wavefront_debug_epilogue(
    wavefront_aligner_t* const wf_aligner) {
  // Print Summary
  if (wf_aligner->system.verbose >= 1) {
    timer_stop(&wf_aligner->system.timer);
    if (wf_aligner->system.verbose >= 4) {
      wavefront_report_verbose_end(stderr,wf_aligner);
    }
    wavefront_report_lite(stderr,wf_aligner);
  }
}
/*
 * Check
 */
void wavefront_debug_check_correct(
    wavefront_aligner_t* const wf_aligner) {
  // Check correct
  if (wf_aligner->system.check_alignment_correct &&
      wf_aligner->align_status.status == WF_STATUS_SUCCESSFUL &&
      wf_aligner->alignment_scope == compute_alignment) {
    if (!wavefront_check_alignment(stderr,wf_aligner)) {
      fprintf(stderr,"[WFA::Check] Error: Alignment incorrect\n");
      exit(1);
    }
  }
}




