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
 * DESCRIPTION: WaveFront penalties handling module
 */

#ifndef WAVEFRONT_WAVEFRONT_PENALTIES_H_
#define WAVEFRONT_WAVEFRONT_PENALTIES_H_

#include "../alignment/linear_penalties.h"
#include "../alignment/affine_penalties.h"
#include "../alignment/affine2p_penalties.h"

/*
 * Distance metrics
 */
typedef enum {
  indel         = 0, // Longest Common Subsequence - LCS
  edit          = 1, // Levenshtein
  gap_linear    = 2, // Needleman-Wunsch
  gap_affine    = 3, // Smith-Waterman-Gotoh
  gap_affine_2p = 4  // Gap-Affine 2-pieces
} distance_metric_t;

/*
 * Wavefront Penalties
 */
typedef struct {
  distance_metric_t distance_metric;  // Alignment metric/distance used
  int match;             // (M <= 0) (Internal variable change to M=0 for WFA)
  int mismatch;          // (X > 0)
  int gap_opening1;      // (O1 >= 0)
  int gap_extension1;    // (E1 > 0)
  int gap_opening2;      // (O2 >= 0)
  int gap_extension2;    // (E2 > 0)
} wavefront_penalties_t;

/*
 * Compute SW-score equivalent (thanks to Eizenga's formula)
 */
#define WF_SCORE_TO_SW_SCORE(swg_match,plen,tlen,wf_score) ((swg_match*(plen+tlen) - wf_score)/2)

/*
 * Penalties adjustment
 */
void wavefront_penalties_set_indel(
    wavefront_penalties_t* const wf_penalties);
void wavefront_penalties_set_edit(
    wavefront_penalties_t* const wf_penalties);
void wavefront_penalties_set_linear(
    wavefront_penalties_t* const wf_penalties,
    linear_penalties_t* const linear_penalties);
void wavefront_penalties_set_affine(
    wavefront_penalties_t* const wf_penalties,
    affine_penalties_t* const affine_penalties);
void wavefront_penalties_set_affine2p(
    wavefront_penalties_t* const wf_penalties,
    affine2p_penalties_t* const affine2p_penalties);

/*
 * Score conversion
 */
int wavefront_penalties_get_score_indel(
    wavefront_penalties_t* const wf_penalties,
    const int score);
int wavefront_penalties_get_score_edit(
    wavefront_penalties_t* const wf_penalties,
    const int score);
int wavefront_penalties_get_score_linear(
    wavefront_penalties_t* const wf_penalties,
    const int score);
int wavefront_penalties_get_score_affine(
    wavefront_penalties_t* const wf_penalties,
    const int score);
int wavefront_penalties_get_score_affine2p(
    wavefront_penalties_t* const wf_penalties,
    const int score);

/*
 * Display
 */
void wavefront_penalties_print(
    FILE* const stream,
    wavefront_penalties_t* const wf_penalties);

#endif /* WAVEFRONT_WAVEFRONT_PENALTIES_H_ */
