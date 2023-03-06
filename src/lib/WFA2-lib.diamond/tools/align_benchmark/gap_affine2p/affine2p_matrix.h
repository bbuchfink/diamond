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

#ifndef AFFINE2P_MATRIX_H_
#define AFFINE2P_MATRIX_H_

#include "utils/commons.h"
#include "alignment/cigar.h"

/*
 * Constants
 */
#define AFFINE2P_SCORE_MAX (10000000)

/*
 * Affine 2-piece Matrix
 */
typedef struct {
  int M;  // Alignment matching/mismatching
  int I1; // Alignment ends with a gap in the reference (insertion)
  int D1; // Alignment ends with a gap in the read (deletion)
  int I2; // Alignment ends with a gap in the reference (insertion)
  int D2; // Alignment ends with a gap in the read (deletion)
} affine2p_cell_t;
typedef struct {
  // Matrix
  affine2p_cell_t** columns;
  int num_rows;
  int num_columns;
} affine2p_matrix_t;

/*
 * Gap-Affine Matrix Setup
 */
void affine2p_matrix_allocate(
    affine2p_matrix_t* const matrix,
    const int num_rows,
    const int num_columns,
    mm_allocator_t* const mm_allocator);
void affine2p_matrix_free(
    affine2p_matrix_t* const matrix,
    mm_allocator_t* const mm_allocator);

/*
 * Display
 */
void affine2p_matrix_print(
    FILE* const stream,
    const affine2p_matrix_t* const matrix,
    const char* const pattern,
    const char* const text);
void affine2p_matrix_print_extended(
    FILE* const stream,
    const affine2p_matrix_t* const matrix,
    const char* const pattern,
    const char* const text);

#endif /* AFFINE2P_MATRIX_H_ */
