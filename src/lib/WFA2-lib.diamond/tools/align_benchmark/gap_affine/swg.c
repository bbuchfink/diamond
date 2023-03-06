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
 * DESCRIPTION: Dynamic-programming algorithm for computing gap-affine
 *   pairwise alignment (Smith-Waterman-Gotoh - SWG)
 */

#include "utils/commons.h"
#include "system/mm_allocator.h"
#include "gap_affine/swg.h"

/*
 * SWG traceback
 */
void swg_traceback(
    affine_matrix_t* const affine_matrix,
    affine_penalties_t* const penalties,
    const int pattern_length,
    const int text_length,
    const int target_v,
    const int target_h,
    cigar_t* const cigar) {
  // Parameters
  affine_cell_t** const dp = affine_matrix->columns;
  char* const operations = cigar->operations;
  cigar->end_offset = cigar->max_operations;
  int op_sentinel = cigar->end_offset - 1;
  // Add final insertions/deletions
  int i;
  for (i=target_v;i<pattern_length;++i) { operations[op_sentinel--] = 'D'; }
  for (i=target_h;i<text_length;++i) { operations[op_sentinel--] = 'I'; }
  // Compute traceback
  affine_matrix_type matrix_type = affine_matrix_M;
  int h = target_h;
  int v = target_v;
  while (h>0 && v>0) {
    switch (matrix_type) {
      case affine_matrix_D:
        // Traceback D-matrix
        operations[op_sentinel--] = 'D';
        if (dp[h][v].D != dp[h][v-1].D + penalties->gap_extension) {
          matrix_type = affine_matrix_M;
        }
        --v;
        break;
      case affine_matrix_I:
        // Traceback I-matrix
        operations[op_sentinel--] = 'I';
        if (dp[h][v].I != dp[h-1][v].I + penalties->gap_extension) {
          matrix_type = affine_matrix_M;
        }
        --h;
        break;
      case affine_matrix_M:
        // Traceback M-matrix
        if (dp[h][v].M == dp[h-1][v-1].M + penalties->mismatch) {
          operations[op_sentinel--] = 'X';
          --h; --v;
        } else if (dp[h][v].M == dp[h][v].D) {
          matrix_type = affine_matrix_D;
        } else if (dp[h][v].M == dp[h][v].I) {
          matrix_type = affine_matrix_I;
        } else if (dp[h][v].M == dp[h-1][v-1].M + penalties->match) {
          operations[op_sentinel--] = 'M';
          --h; --v;
        } else {
          fprintf(stderr,"SWG backtrace. No backtrace operation found\n");
          exit(1);
        }
        break;
    }
  }
  // Add initial deletions/insertions
  while (v>0) { operations[op_sentinel--] = 'D'; --v; }
  while (h>0) { operations[op_sentinel--] = 'I'; --h; }
  cigar->begin_offset = op_sentinel+1;
  cigar->score = dp[target_h][target_v].M;
}
/*
 * SWG alignment
 */
void swg_align(
    affine_matrix_t* const affine_matrix,
    affine_penalties_t* const penalties,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar) {
  // Parameters
  affine_cell_t** const dp = affine_matrix->columns;
  int h, v;
  // Init DP
  dp[0][0].D = AFFINE_SCORE_MAX;
  dp[0][0].I = AFFINE_SCORE_MAX;
  dp[0][0].M = 0;
  for (v=1;v<=pattern_length;++v) { // Init first column
    dp[0][v].D = penalties->gap_opening + v*penalties->gap_extension;
    dp[0][v].I = AFFINE_SCORE_MAX;
    dp[0][v].M = dp[0][v].D;
  }
  for (h=1;h<=text_length;++h) { // Init first row
    dp[h][0].D = AFFINE_SCORE_MAX;
    dp[h][0].I = penalties->gap_opening + h*penalties->gap_extension;
    dp[h][0].M = dp[h][0].I;
  }
  // Compute DP
  for (h=1;h<=text_length;++h) {
    for (v=1;v<=pattern_length;++v) {
      // Update DP.D
      const int del_new = dp[h][v-1].M + penalties->gap_opening + penalties->gap_extension;
      const int del_ext = dp[h][v-1].D + penalties->gap_extension;
      const int del = MIN(del_new,del_ext);
      dp[h][v].D = del;
      // Update DP.I
      const int ins_new = dp[h-1][v].M + penalties->gap_opening + penalties->gap_extension;
      const int ins_ext = dp[h-1][v].I + penalties->gap_extension;
      const int ins = MIN(ins_new,ins_ext);
      dp[h][v].I = ins;
      // Update DP.M
      const int m_match = dp[h-1][v-1].M + ((pattern[v-1]==text[h-1]) ? penalties->match : penalties->mismatch);
      dp[h][v].M = MIN(m_match,MIN(ins,del));
    }
  }
  // Compute traceback
  swg_traceback(
      affine_matrix,penalties,
      pattern_length,text_length,
      pattern_length,text_length,cigar);
  // DEBUG
  //affine_matrix_print(stderr,affine_matrix,pattern,text);
}
/*
 * SWG alignment (ends-free)
 */
void swg_align_endsfree(
    affine_matrix_t* const affine_matrix,
    affine_penalties_t* const penalties,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free,
    cigar_t* const cigar) {
  // Parameters
  affine_cell_t** const dp = affine_matrix->columns;
  const int pattern_min_v = pattern_length - pattern_end_free;
  const int text_min_h = text_length - text_end_free;
  int h, v;
  // Init DP
  dp[0][0].D = AFFINE_SCORE_MAX;
  dp[0][0].I = AFFINE_SCORE_MAX;
  dp[0][0].M = 0;
  for (v=1;v<=pattern_length;++v) { // Init first column
    dp[0][v].D = (v > pattern_begin_free) ? penalties->gap_opening + (v+1-pattern_begin_free)*penalties->gap_extension : 0;
    dp[0][v].I = AFFINE_SCORE_MAX;
    dp[0][v].M = dp[0][v].D;
  }
  for (h=1;h<=text_length;++h) { // Init first row
    dp[h][0].D = AFFINE_SCORE_MAX;
    dp[h][0].I = (h > text_begin_free) ? penalties->gap_opening + (h+1-text_begin_free)*penalties->gap_extension : 0;
    dp[h][0].M = dp[h][0].I;
  }
  // Keep minimum score
  int min_v = 0, min_h = 0, min_score = AFFINE_SCORE_MAX;
  if (text_min_h==0 && dp[0][pattern_length].M<min_score) {
    min_score = dp[0][pattern_length].M;
    min_v = pattern_length;
    min_h = 0;
  }
  // Compute DP
  for (h=1;h<=text_length;++h) {
    for (v=1;v<=pattern_length;++v) {
      // Update DP.D
      const int del_new = dp[h][v-1].M + penalties->gap_opening + penalties->gap_extension;
      const int del_ext = dp[h][v-1].D + penalties->gap_extension;
      const int del = MIN(del_new,del_ext);
      dp[h][v].D = del;
      // Update DP.I
      const int ins_new = dp[h-1][v].M + penalties->gap_opening + penalties->gap_extension;
      const int ins_ext = dp[h-1][v].I + penalties->gap_extension;
      const int ins = MIN(ins_new,ins_ext);
      dp[h][v].I = ins;
      // Update DP.M
      const int m_match = dp[h-1][v-1].M + ((pattern[v-1]==text[h-1]) ? penalties->match : penalties->mismatch);
      dp[h][v].M = MIN(m_match,MIN(ins,del));
    }
    // Pattern aligned. Keep minimum score
    if (h>=text_min_h && dp[h][pattern_length].M<min_score) {
      min_score = dp[h][pattern_length].M;
      min_v = pattern_length;
      min_h = h;
    }
  }
  // Text aligned. Keep minimum score
  for (v=pattern_length;v>=0;--v) {
    if (v>=pattern_min_v && dp[text_length][v].M<min_score) {
      min_score = dp[text_length][v].M;
      min_v = v;
      min_h = text_length;
    }
  }
  // Compute traceback
  swg_traceback(
      affine_matrix,penalties,
      pattern_length,text_length,
      min_v,min_h,cigar);
  // DEBUG
  //affine_matrix_print(stderr,affine_matrix,pattern,text);
}
/*
 * SWG alignment (banded)
 */
void swg_align_banded(
    affine_matrix_t* const affine_matrix,
    affine_penalties_t* const penalties,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int bandwidth,
    cigar_t* const cigar) {
  // Parameters
  const int k_end = ABS(text_length-pattern_length)+1;
  int effective_bandwidth = bandwidth;
  if (effective_bandwidth < k_end) effective_bandwidth = k_end;
  if (effective_bandwidth > pattern_length) effective_bandwidth = pattern_length;
  if (effective_bandwidth > text_length) effective_bandwidth = text_length;
  affine_cell_t** const dp = affine_matrix->columns;
  int h, v;
  // Initialize
  dp[0][0].D = AFFINE_SCORE_MAX;
  dp[0][0].I = AFFINE_SCORE_MAX;
  dp[0][0].M = 0;
  for (v=1;v<=effective_bandwidth;++v) {
    dp[0][v].D = penalties->gap_opening + v*penalties->gap_extension;
    dp[0][v].I = AFFINE_SCORE_MAX;
    dp[0][v].M = dp[0][v].D;
  }
  // Compute DP
  for (h=1;h<=text_length;++h) {
    // Compute lo limit
    int lo;
    if (h <= effective_bandwidth) {
      lo = 1;
      dp[h][lo-1].D = AFFINE_SCORE_MAX;
      dp[h][lo-1].I = penalties->gap_opening + h*penalties->gap_extension;
      dp[h][lo-1].M = dp[h][lo-1].I;
    } else {
      lo = h - effective_bandwidth;
      dp[h][lo-1].D = AFFINE_SCORE_MAX;
      dp[h][lo-1].I = AFFINE_SCORE_MAX;
      dp[h][lo-1].M = AFFINE_SCORE_MAX;
    }
    // Compute hi limit
    int hi = h + effective_bandwidth - 1;
    if (hi > pattern_length) {
      hi = pattern_length;
    } else if (h > 1) {
      dp[h-1][hi].D = AFFINE_SCORE_MAX;
      dp[h-1][hi].I = AFFINE_SCORE_MAX;
      dp[h-1][hi].M = AFFINE_SCORE_MAX;
    }
    // Compute column
    for (v=lo;v<=hi;++v) {
      // Update DP.D
      const int del_new = dp[h][v-1].M + penalties->gap_opening + penalties->gap_extension;
      const int del_ext = dp[h][v-1].D + penalties->gap_extension;
      const int del = MIN(del_new,del_ext);
      dp[h][v].D = del;
      // Update DP.I
      const int ins_new = dp[h-1][v].M + penalties->gap_opening + penalties->gap_extension;
      const int ins_ext = dp[h-1][v].I + penalties->gap_extension;
      const int ins = MIN(ins_new,ins_ext);
      dp[h][v].I = ins;
      // Update DP.M
      const int m_match = dp[h-1][v-1].M + ((pattern[v-1]==text[h-1]) ? penalties->match : penalties->mismatch);
      dp[h][v].M = MIN(m_match,MIN(ins,del));
    }
  }
  // Compute traceback
  swg_traceback(
      affine_matrix,penalties,
      pattern_length,text_length,
      pattern_length,text_length,cigar);
}
