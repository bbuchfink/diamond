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
 * DESCRIPTION: Dynamic-programming algorithm to compute indel alignment (LCS)
 */

#include "indel/indel_dp.h"

/*
 * Indel distance computation using dynamic-programming matrix
 */
void indel_dp_traceback(
    score_matrix_t* const score_matrix,
    cigar_t* const cigar) {
  // Parameters
  int** const dp = score_matrix->columns;
  char* const operations = cigar->operations;
  cigar->end_offset = cigar->max_operations;
  int op_sentinel = cigar->end_offset - 1;
  int h, v;
  // Compute traceback
  h = score_matrix->num_columns-1;
  v = score_matrix->num_rows-1;
  while (h>0 && v>0) {
    if (dp[h][v]==dp[h][v-1]+1) {
      operations[op_sentinel--] = 'D';
      --v;
    } else if (dp[h][v]==dp[h-1][v]+1) {
      operations[op_sentinel--] = 'I';
      --h;
    } else if (dp[h][v]==dp[h-1][v-1]) {
      operations[op_sentinel--] = 'M';
      --h;
      --v;
    } else {
      fprintf(stderr,"Indel backtrace. No backtrace operation found");
      exit(1);
    }
  }
  while (h>0) {operations[op_sentinel--] = 'I'; --h;}
  while (v>0) {operations[op_sentinel--] = 'D'; --v;}
  cigar->begin_offset = op_sentinel + 1;
}
void indel_dp_compute(
    score_matrix_t* const score_matrix,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar) {
  // Parameters
  int** dp = score_matrix->columns;
  int h, v;
  // Init DP
  for (v=0;v<=pattern_length;++v) dp[0][v] = v; // No ends-free
  for (h=0;h<=text_length;++h) dp[h][0] = h; // No ends-free
  // Compute DP
  for (h=1;h<=text_length;++h) {
    for (v=1;v<=pattern_length;++v) {
      int min = MIN(dp[h-1][v]+1,dp[h][v-1]+1); // Indel
      if (text[h-1] == pattern[v-1]) {
        min = MIN(min,dp[h-1][v-1]);
      }
      dp[h][v] = min;
    }
  }
  // Compute traceback
  indel_dp_traceback(score_matrix,cigar);
  // DEBUG
  // score_matrix_print(stderr,score_matrix,pattern,text);
}
