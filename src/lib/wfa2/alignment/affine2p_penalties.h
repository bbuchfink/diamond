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
 * DESCRIPTION: Gap-Affine 2-Pieces penalties
 */

#ifndef AFFINE2P_PENALTIES_H_
#define AFFINE2P_PENALTIES_H_

#include "../utils/commons.h"

/*
 * Affine 2-piece penalties
 */
typedef struct {
  int match;             // (Penalty representation; usually M <= 0)
  int mismatch;          // (Penalty representation; usually X > 0)
  // Usually concave; Q1 + E1 < Q2 + E2 and E1 > E2.
  int gap_opening1;      // (Penalty representation; usually O1 > 0)
  int gap_extension1;    // (Penalty representation; usually E1 > 0)
  int gap_opening2;      // (Penalty representation; usually O2 > 0)
  int gap_extension2;    // (Penalty representation; usually E2 > 0)
} affine2p_penalties_t;

/*
 * Affine 2-piece matrix-type (for bcktrace)
 */
typedef enum {
  affine2p_matrix_M,
  affine2p_matrix_I1,
  affine2p_matrix_I2,
  affine2p_matrix_D1,
  affine2p_matrix_D2
} affine2p_matrix_type;

#endif /* AFFINE2P_PENALTIES_H_ */
