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
 * DESCRIPTION: Support functions for wavefront heuristic strategies
 */

#ifndef WAVEFRONT_HEURISTIC_H_
#define WAVEFRONT_HEURISTIC_H_

#include "../utils/commons.h"

// Wavefront ahead definition
typedef struct _wavefront_aligner_t wavefront_aligner_t;

/*
 * Wavefront Heuristics
 */
typedef enum {
  wf_heuristic_none            = 0x0000000000000000ul,
  wf_heuristic_banded_static   = 0x0000000000000001ul,
  wf_heuristic_banded_adaptive = 0x0000000000000002ul,
  wf_heuristic_wfadaptive      = 0x0000000000000004ul,
  wf_heuristic_xdrop           = 0x0000000000000010ul,
  wf_heuristic_zdrop           = 0x0000000000000020ul,
  wf_heuristic_wfmash          = 0x0000000000000040ul,
} wf_heuristic_strategy;
typedef struct {
  // Heuristic
  wf_heuristic_strategy strategy;     // Heuristic strategy
  int steps_between_cutoffs;          // Score-steps between heuristic cut-offs
  // Static/Adaptive Banded
  int min_k;                          // Banded: Minimum k to consider in band
  int max_k;                          // Banded: Maximum k to consider in band
  // WFAdaptive
  int min_wavefront_length;           // Adaptive: Minimum wavefronts length to cut-off
  int max_distance_threshold;         // Adaptive: Maximum distance between offsets allowed
  // Drops
  int xdrop;                          // X-drop parameter
  int zdrop;                          // Z-drop parameter
  // Internals
  int steps_wait;                     // Score-steps until next cut-off
  int max_sw_score;                   // Maximum score observed (for x/z drops)
  int max_sw_score_offset;            // Offset of the maximum score observed
  int max_sw_score_k;                 // Diagonal of the maximum score observed
} wavefront_heuristic_t;

/*
 * Setup
 */
void wavefront_heuristic_set_none(
    wavefront_heuristic_t* const wf_heuristic);

void wavefront_heuristic_set_wfadaptive(
    wavefront_heuristic_t* const wf_heuristic,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int steps_between_cutoffs);
void wavefront_heuristic_set_wfmash(
    wavefront_heuristic_t* const wf_heuristic,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int steps_between_cutoffs);

void wavefront_heuristic_set_xdrop(
    wavefront_heuristic_t* const wf_heuristic,
    const int xdrop,
    const int steps_between_cutoffs);
void wavefront_heuristic_set_zdrop(
    wavefront_heuristic_t* const wf_heuristic,
    const int ydrop,
    const int steps_between_cutoffs);

void wavefront_heuristic_set_banded_static(
    wavefront_heuristic_t* const wf_heuristic,
    const int band_min_k,
    const int band_max_k);
void wavefront_heuristic_set_banded_adaptive(
    wavefront_heuristic_t* const wf_heuristic,
    const int band_min_k,
    const int band_max_k,
    const int steps_between_cutoffs);

void wavefront_heuristic_clear(
    wavefront_heuristic_t* const wf_heuristic);

/*
 * Wavefront heuristic cut-off
 */
void wavefront_heuristic_cufoff(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int score_mod);

/*
 * Display
 */
void wavefront_heuristic_print(
    FILE* const stream,
    wavefront_heuristic_t* const wf_heuristic);

#endif /* WAVEFRONT_HEURISTIC_H_ */
