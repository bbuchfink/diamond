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
 * DESCRIPTION: WaveFront alignment module for computing wavefronts (gap-affine)
 */

#ifndef WAVEFRONT_COMPUTE_AFFINE_H_
#define WAVEFRONT_COMPUTE_AFFINE_H_

#include "wavefront_aligner.h"

/*
 * Compute Kernels
 */
void wavefront_compute_affine_idm(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi);
void wavefront_compute_affine_idm_piggyback(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi);

/*
 * Compute Wavefronts (gap-affine)
 */
void wavefront_compute_affine(
    wavefront_aligner_t* const wf_aligner,
    const int score);

#endif /* WAVEFRONT_COMPUTE_AFFINE_H_ */
