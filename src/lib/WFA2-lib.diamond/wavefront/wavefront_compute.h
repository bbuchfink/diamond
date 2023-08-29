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
 * DESCRIPTION: WaveFront alignment module for computing wavefronts
 */

#ifndef WAVEFRONT_COMPUTE_H_
#define WAVEFRONT_COMPUTE_H_

#include "wavefront_aligner.h"

/*
 * Compute limits
 */
void wavefront_compute_limits_input(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    int* const lo,
    int* const hi);
void wavefront_compute_limits_output(
    wavefront_aligner_t* const wf_aligner,
    const int lo,
    const int hi,
    int* const effective_lo,
    int* const effective_hi);

/*
 * Score translation
 */
int wavefront_compute_classic_score(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length,
    const int wf_score);

/*
 * Input wavefronts (fetch)
 */
void wavefront_compute_fetch_input(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set,
    const int score);

/*
 * Output wavefronts (allocate)
 */
void wavefront_compute_allocate_output_null(
    wavefront_aligner_t* const wf_aligner,
    const int score);
void wavefront_compute_allocate_output(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set,
    const int score,
    const int lo,
    const int hi);

/*
 * Initialize wavefronts ends
 */
void wavefront_compute_init_ends(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi);

/*
 * Process wavefronts ends
 */
void wavefront_compute_trim_ends(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront);
void wavefront_compute_process_ends(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set,
    const int score);

/*
 * Multithread dispatcher
 */
#ifdef WFA_PARALLEL
int wavefront_compute_num_threads(
    wavefront_aligner_t* const wf_aligner,
    const int lo,
    const int hi);
void wavefront_compute_thread_limits(
    const int thread_id,
    const int num_theads,
    const int lo,
    const int hi,
    int* const thread_lo,
    int* const thread_hi);
#else
#define wavefront_compute_num_threads(wf_aligner,lo,hi) 1
#endif

#endif /* WAVEFRONT_COMPUTE_H_ */
