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
 * DESCRIPTION: Gap-affine Matrix
 */

#ifndef AFFINE_MATRIX_H_
#define AFFINE_MATRIX_H_

#include "alignment/cigar.h"

/*
 * Constants
 */
#define AFFINE_SCORE_MAX (10000000)

/*
 * Affine Matrix
 */
typedef struct {
  int M; // Alignment matching/mismatching
  int I; // Alignment ends with a gap in the reference (insertion)
  int D; // Alignment ends with a gap in the read (deletion)
} affine_cell_t;
typedef struct {
  // Affine Matrix
  affine_cell_t** columns;
  int num_rows;
  int num_columns;
} affine_matrix_t;

/*
 * Gap-Affine Matrix Setup
 */
void affine_matrix_allocate(
    affine_matrix_t* const matrix,
    const int num_rows,
    const int num_columns,
    mm_allocator_t* const mm_allocator);
void affine_matrix_free(
    affine_matrix_t* const matrix,
    mm_allocator_t* const mm_allocator);

/*
 * Display
 */
void affine_matrix_print(
    FILE* const stream,
    const affine_matrix_t* const matrix,
    const char* const pattern,
    const char* const text);
void affine_matrix_print_extended(
    FILE* const stream,
    const affine_matrix_t* const matrix,
    const char* const pattern,
    const char* const text);

#endif /* AFFINE_MATRIX_H_ */
