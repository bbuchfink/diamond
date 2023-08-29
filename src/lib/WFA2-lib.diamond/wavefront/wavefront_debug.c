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

#include "utils/commons.h"
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
  wavefront_sequences_t* const sequences = (wf_aligner->bialigner==NULL) ?
      &wf_aligner->sequences : &wf_aligner->bialigner->alg_forward->sequences;
  const char* const pattern = sequences->pattern_buffer;
  const int pattern_length = sequences->pattern_buffer_length;
  const char* const text = sequences->text_buffer;
  const int text_length = sequences->text_buffer_length;
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
        if (sequences->mode != wf_sequences_lambda) {
          const bool is_match = (pattern[pattern_pos]==text[text_pos]);
          if (!is_match) {
            fprintf(stream,"[WFA::Check] Alignment not matching (pattern[%d]=%c != text[%d]=%c)\n",
                pattern_pos,pattern[pattern_pos],
                text_pos,text[text_pos]);
            alignment_correct = false;
            break;
          }
        }
        ++pattern_pos;
        ++text_pos;
        break;
      }
      case 'X': {
        // Check mismatch
        if (sequences->mode != wf_sequences_lambda) {
          const bool is_match = (pattern[pattern_pos]==text[text_pos]);
          if (is_match) {
            fprintf(stream,"[WFA::Check] Alignment not mismatching (pattern[%d]=%c == text[%d]=%c)\n",
                pattern_pos,pattern[pattern_pos],
                text_pos,text[text_pos]);
            alignment_correct = false;
            break;
          }
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
  wavefront_sequences_t* const sequences = (wf_aligner->bialigner==NULL) ?
      &wf_aligner->sequences : &wf_aligner->bialigner->alg_subsidiary->sequences;
  const char* const pattern = sequences->pattern;
  const int pattern_length = sequences->pattern_length;
  const char* const text = sequences->text;
  const int text_length = sequences->text_length;
  const int status = wf_aligner->align_status.status;
  const uint64_t memory_used = wf_aligner->align_status.memory_used;
  // BANNER (#0)
  fprintf(stream,"[WFA::Debug]");
  // SCORE (#1)
  //  const int score = wavefront_compute_classic_score(
  //      wf_aligner,pattern_length,text_length,wf_aligner->cigar->score);
  const int score = wf_aligner->cigar->score;
  fprintf(stream,"\t%d",(score==INT32_MIN) ? -1 : score);
  // PATTERN_LENGTH (#2)
  fprintf(stream,"\t%d",pattern_length);
  // TEXT_LENGTH (#3)
  fprintf(stream,"\t%d",text_length);
  // STATUS (#4)
  fprintf(stream,"\t%s",wavefront_align_strerror_short(status));
  // TIME (#5)
  fprintf(stream,"\t%2.3f",TIMER_GET_TOTAL_MS(&wf_aligner->system.timer));
  // MEMORY (#6)
  fprintf(stream,"\t%luMB\t",CONVERT_B_TO_MB(memory_used));
  // ATTRIBUTES (#7)
  fprintf(stream,"[");
  fprintf(stream,"%d",wf_aligner->align_status.status);
  fprintf(stream,";");
  wavefront_aligner_print_mode(stream,wf_aligner);
  fprintf(stream,";");
  wavefront_aligner_print_scope(stream,wf_aligner);
  fprintf(stream,";");
  wavefront_penalties_print(stream,&wf_aligner->penalties);
  fprintf(stream,";");
  wavefront_aligner_print_conf(stream,wf_aligner);
  fprintf(stream,";");
  wavefront_heuristic_print(stream,&wf_aligner->heuristic);
  fprintf(stream,";");
  fprintf(stream,"(%d,%d,%d)",
      wf_aligner->wf_components.num_wavefronts,
      wf_aligner->wf_components.historic_min_lo,
      wf_aligner->wf_components.historic_max_hi);
  fprintf(stream,"]\t");
  // CIGAR (#8)
  if (cigar_is_null(wf_aligner->cigar)) {
    fprintf(stream,"-");
  } else {
    cigar_print(stream,wf_aligner->cigar,true);
  }
  // SEQUENCES (#9 #10)
  if (sequences->mode == wf_sequences_lambda) {
    fprintf(stream,"\t-\t-");
  } else {
    fprintf(stream,"\t%.*s\t%.*s",pattern_length,pattern,text_length,text);
  }
  fprintf(stream,"\n");
}
/*
 * Debug
 */
void wavefront_debug_begin(
    wavefront_aligner_t* const wf_aligner) {
  // Check verbose level
  if (wf_aligner->system.verbose >= 1) {
    timer_reset(&wf_aligner->system.timer);
    timer_start(&wf_aligner->system.timer);
  }
}
void wavefront_debug_end(
    wavefront_aligner_t* const wf_aligner) {
  // Print Summary
  if (wf_aligner->system.verbose >= 1) {
    timer_stop(&wf_aligner->system.timer);
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




