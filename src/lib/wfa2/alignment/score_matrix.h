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
 * DESCRIPTION: Score matrix for alignment using dynamic programming
 */

#ifndef SCORE_MATRIX_H_
#define SCORE_MATRIX_H_

#include "../utils/commons.h"
#include "../system/mm_allocator.h"
#include "../alignment/cigar.h"

/*
 * Constants
 */
#define SCORE_MAX (10000000)

/*
 * Score Matrix
 */
typedef struct {
  // Score Columns
  int** columns;
  int num_rows;
  int num_columns;
  // MM
  mm_allocator_t* mm_allocator;
} score_matrix_t;

/*
 * Setup
 */
void score_matrix_allocate(
    score_matrix_t* const score_matrix,
    const int num_rows,
    const int num_columns,
    mm_allocator_t* const mm_allocator);
void score_matrix_free(
    score_matrix_t* const score_matrix);

/*
 * Display
 */
void score_matrix_print(
    FILE* const stream,
    const score_matrix_t* const score_matrix,
    const char* const pattern,
    const char* const text);

#endif /* SCORE_MATRIX_H_ */
