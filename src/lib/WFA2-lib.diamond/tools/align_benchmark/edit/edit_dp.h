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
 * DESCRIPTION: Dynamic-programming algorithm to compute Levenshtein alignment (edit)
 */

#ifndef EDIT_DP_H_
#define EDIT_DP_H_

#include "alignment/score_matrix.h"

/*
 * Edit distance computation using dynamic-programming matrix
 */
void edit_dp_align(
    score_matrix_t* const score_matrix,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar);
/*
 * Edit distance computation using dynamic-programming matrix (banded)
 */
void edit_dp_align_banded(
    score_matrix_t* const score_matrix,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int bandwidth,
    cigar_t* const cigar);

#endif /* EDIT_DP_H_ */
