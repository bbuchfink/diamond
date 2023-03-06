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
 * DESCRIPTION: Dynamic-programming algorithm for computing
 *   pairwise alignment using gap-affine 2-piece penalties
 */

#include "affine2p_dp.h"

/*
 * Gap-affine 2-piece backtrace using dynamic-programming matrix
 */
void affine2p_dp_traceback(
    affine2p_matrix_t* const matrix,
    affine2p_penalties_t* const penalties,
    const int pattern_length,
    const int text_length,
    cigar_t* const cigar) {
  // Parameters
  affine2p_cell_t** const dp = matrix->columns;
  char* const operations = cigar->operations;
  cigar->end_offset = cigar->max_operations;
  int op_sentinel = cigar->end_offset - 1;
  // Compute traceback
  int h = text_length;
  int v = pattern_length;
  affine2p_matrix_type matrix_type = affine2p_matrix_M;
  while (h>0 && v>0) {
    switch (matrix_type) {
      case affine2p_matrix_D1:
        // Traceback D1-matrix
        operations[op_sentinel--] = 'D';
        if (dp[h][v].D1 != dp[h][v-1].D1 + penalties->gap_extension1) {
          matrix_type = affine2p_matrix_M;
        }
        --v;
        break;
      case affine2p_matrix_D2:
        // Traceback D2-matrix
        operations[op_sentinel--] = 'D';
        if (dp[h][v].D2 != dp[h][v-1].D2 + penalties->gap_extension2) {
          matrix_type = affine2p_matrix_M;
        }
        --v;
        break;
      case affine2p_matrix_I1:
        // Traceback I1-matrix
        operations[op_sentinel--] = 'I';
        if (dp[h][v].I1 != dp[h-1][v].I1 + penalties->gap_extension1) {
          matrix_type = affine2p_matrix_M;
        }
        --h;
        break;
      case affine2p_matrix_I2:
        // Traceback I2-matrix
        operations[op_sentinel--] = 'I';
        if (dp[h][v].I2 != dp[h-1][v].I2 + penalties->gap_extension2) {
          matrix_type = affine2p_matrix_M;
        }
        --h;
        break;
      case affine2p_matrix_M:
        // Traceback M-matrix
        if (dp[h][v].M == dp[h-1][v-1].M + penalties->mismatch) {
          operations[op_sentinel--] = 'X';
          --h; --v;
        } else if (dp[h][v].M == dp[h][v].D2) {
          matrix_type = affine2p_matrix_D2;
        } else if (dp[h][v].M == dp[h][v].D1) {
          matrix_type = affine2p_matrix_D1;
        } else if (dp[h][v].M == dp[h][v].I2) {
          matrix_type = affine2p_matrix_I2;
        } else if (dp[h][v].M == dp[h][v].I1) {
          matrix_type = affine2p_matrix_I1;
        } else if (dp[h][v].M == dp[h-1][v-1].M + penalties->match) {
          operations[op_sentinel--] = 'M';
          --h; --v;
        } else {
          fprintf(stderr,"Affine-2P Backtrace. No backtrace operation found\n");
          exit(1);
        }
        break;
    }
  }
  while (h>0) {operations[op_sentinel--] = 'I'; --h;}
  while (v>0) {operations[op_sentinel--] = 'D'; --v;}
  cigar->begin_offset = op_sentinel+1;
}
/*
 * Gap-affine 2-piece alignment computation using dynamic-programming matrix
 */
void affine2p_dp_align(
    affine2p_matrix_t* const matrix,
    affine2p_penalties_t* const penalties,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar) {
  // Parameters
  affine2p_cell_t** const dp = matrix->columns;
  int h, v;
  // Initial conditions
  dp[0][0].I1 = AFFINE2P_SCORE_MAX;
  dp[0][0].I2 = AFFINE2P_SCORE_MAX;
  dp[0][0].D1 = AFFINE2P_SCORE_MAX;
  dp[0][0].D2 = AFFINE2P_SCORE_MAX;
  dp[0][0].M = 0;
  for (v=1;v<=pattern_length;++v) { // Initialize first column
    dp[0][v].I1 = AFFINE2P_SCORE_MAX;
    dp[0][v].I2 = AFFINE2P_SCORE_MAX;
    dp[0][v].D1 = penalties->gap_opening1 + v*penalties->gap_extension1;
    dp[0][v].D2 = penalties->gap_opening2 + v*penalties->gap_extension2;
    dp[0][v].M = MIN(dp[0][v].D1,dp[0][v].D2);
  }
  for (h=1;h<=text_length;++h) { // Initialize first row
    dp[h][0].I1 = penalties->gap_opening1 + h*penalties->gap_extension1;
    dp[h][0].I2 = penalties->gap_opening2 + h*penalties->gap_extension2;
    dp[h][0].D1 = AFFINE2P_SCORE_MAX;
    dp[h][0].D2 = AFFINE2P_SCORE_MAX;
    dp[h][0].M = MIN(dp[h][0].I1,dp[h][0].I2);
  }
  // Compute DP
  for (h=1;h<=text_length;++h) {
    for (v=1;v<=pattern_length;++v) {
      // Update DP.I
      const int ins1_new = dp[h-1][v].M + penalties->gap_opening1 + penalties->gap_extension1;
      const int ins1_ext = dp[h-1][v].I1 + penalties->gap_extension1;
      const int ins1 = MIN(ins1_new,ins1_ext);
      dp[h][v].I1 = ins1;
      const int ins2_new = dp[h-1][v].M + penalties->gap_opening2 + penalties->gap_extension2;
      const int ins2_ext = dp[h-1][v].I2 + penalties->gap_extension2;
      const int ins2 = MIN(ins2_new,ins2_ext);
      dp[h][v].I2 = ins2;
      // Update DP.D
      const int del1_new = dp[h][v-1].M + penalties->gap_opening1 + penalties->gap_extension1;
      const int del1_ext = dp[h][v-1].D1 + penalties->gap_extension1;
      const int del1 = MIN(del1_new,del1_ext);
      dp[h][v].D1 = del1;
      const int del2_new = dp[h][v-1].M + penalties->gap_opening2 + penalties->gap_extension2;
      const int del2_ext = dp[h][v-1].D2 + penalties->gap_extension2;
      const int del2 = MIN(del2_new,del2_ext);
      dp[h][v].D2 = del2;
      // Update DP.M
      const int m_match = dp[h-1][v-1].M + ((pattern[v-1]==text[h-1]) ? penalties->match : penalties->mismatch);
      const int min_dels = MIN(del1,del2);
      const int min_inss = MIN(ins1,ins2);
      const int min_gaps = MIN(min_dels,min_inss);
      dp[h][v].M = MIN(m_match,min_gaps);
    }
  }
  // Compute traceback
  affine2p_dp_traceback(matrix,penalties,pattern_length,text_length,cigar);
//  // DEBUG
//  affine2p_matrix_print(stderr,matrix,pattern,text);
}

