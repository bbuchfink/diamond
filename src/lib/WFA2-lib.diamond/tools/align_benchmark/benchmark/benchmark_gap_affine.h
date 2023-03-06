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
 * DESCRIPTION: Gap-affine alignment algorithms wrapper (including WFA)
 */

#ifndef BENCHMARK_GAP_AFFINE_H_
#define BENCHMARK_GAP_AFFINE_H_

#include "wavefront/wavefront_align.h"
#include "benchmark/benchmark_utils.h"

/*
 * Benchmark SWG
 */
void benchmark_gap_affine_swg(
    align_input_t* const align_input,
    affine_penalties_t* const penalties);
void benchmark_gap_affine_swg_banded(
    align_input_t* const align_input,
    affine_penalties_t* const penalties,
    const int bandwidth);
void benchmark_gap_affine_swg_endsfree(
    align_input_t* const align_input,
    affine_penalties_t* const penalties);
void benchmark_gap_affine_wavefront(
    align_input_t* const align_input,
    affine_penalties_t* const penalties);

#endif /* BENCHMARK_GAP_AFFINE_H_ */
