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
 * DESCRIPTION: Gap-Affine 2-Pieces Matrix (for dynamic-programming methods)
 */

#include "affine2p_matrix.h"

/*
 * Gap-Affine Matrix Setup
 */
void affine2p_matrix_allocate(
    affine2p_matrix_t* const matrix,
    const int num_rows,
    const int num_columns,
    mm_allocator_t* const mm_allocator) {
  // Allocate Matrix
  int h;
  matrix->num_rows = num_rows;
  matrix->num_columns = num_columns;
  matrix->columns = mm_allocator_malloc(mm_allocator,num_columns*sizeof(affine2p_cell_t*)); // Columns
  for (h=0;h<num_columns;++h) {
    matrix->columns[h] = mm_allocator_calloc(mm_allocator,num_rows,affine2p_cell_t,false); // Rows
  }
}
void affine2p_matrix_free(
    affine2p_matrix_t* const matrix,
    mm_allocator_t* const mm_allocator) {
  const int num_columns = matrix->num_columns;
  int h;
  for (h=0;h<num_columns;++h) {
    mm_allocator_free(mm_allocator,matrix->columns[h]);
  }
  mm_allocator_free(mm_allocator,matrix->columns);
}
/*
 * Display
 */
void affine2p_matrix_print_value(
    FILE* const stream,
    const int value) {
  if (value >= 0 && value < AFFINE2P_SCORE_MAX) {
    fprintf(stream,"%2d",value);
  } else {
    fprintf(stream," *");
  }
}
void affine2p_matrix_print_cell(
    FILE* const stream,
    affine2p_cell_t* const affine2p_cell) {
  fprintf(stream," ");
  affine2p_matrix_print_value(stream,affine2p_cell->M);
  fprintf(stream,"{");
  affine2p_matrix_print_value(stream,affine2p_cell->I1);
  fprintf(stream,"/");
  affine2p_matrix_print_value(stream,affine2p_cell->I2);
  fprintf(stream,",");
  affine2p_matrix_print_value(stream,affine2p_cell->D1);
  fprintf(stream,"/");
  affine2p_matrix_print_value(stream,affine2p_cell->D2);
  fprintf(stream,"} ");
}
void affine2p_matrix_print(
    FILE* const stream,
    const affine2p_matrix_t* const matrix,
    const char* const pattern,
    const char* const text) {
  // Parameters
  affine2p_cell_t** const dp = matrix->columns;
  int i, j;
  // Print Header
  fprintf(stream,"     ");
  for (i=0;i<matrix->num_columns-1;++i) {
    fprintf(stream," %c ",text[i]);
  }
  fprintf(stream,"\n ");
  for (i=0;i<matrix->num_columns;++i) {
    fprintf(stream," ");
    affine2p_matrix_print_value(stream,dp[i][0].M);
  }
  fprintf(stream,"\n");
  // Print Rows
  for (i=1;i<matrix->num_rows;++i) {
    fprintf(stream,"%c",pattern[i-1]);
    for (j=0;j<matrix->num_columns;++j) {
      fprintf(stream," ");
      affine2p_matrix_print_value(stream,dp[j][i].M);
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}
void affine2p_matrix_print_extended(
    FILE* const stream,
    const affine2p_matrix_t* const matrix,
    const char* const pattern,
    const char* const text) {
  // Parameters
  affine2p_cell_t** const dp = matrix->columns;
  int i, j;
  // Print Header
  fprintf(stream,"         ");
  for (i=0;i<matrix->num_columns-1;++i) {
    fprintf(stream,"     %c     ",text[i]);
  }
  fprintf(stream,"\n ");
  for (i=0;i<matrix->num_columns;++i) {
    affine2p_matrix_print_cell(stderr,&(dp[i][0]));
  }
  fprintf(stream,"\n");
  // Print Rows
  for (i=1;i<matrix->num_rows;++i) {
    fprintf(stream,"%c",pattern[i-1]);
    for (j=0;j<matrix->num_columns;++j) {
      affine2p_matrix_print_cell(stderr,&(dp[j][i]));
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}
