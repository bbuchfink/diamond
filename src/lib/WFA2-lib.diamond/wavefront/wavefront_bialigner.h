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
 */

#ifndef WAVEFRONT_BIALIGNER_H_
#define WAVEFRONT_BIALIGNER_H_

#include "wavefront_penalties.h"
#include "wavefront_attributes.h"
#include "wavefront_heuristic.h"
#include "wavefront_offset.h"
#include "wavefront_sequences.h"

// Wavefront ahead definition
typedef struct _wavefront_aligner_t wavefront_aligner_t;

typedef struct {
  // Scores
  int score;                      // Score total
  int score_forward;              // Score (forward)
  int score_reverse;              // Score (reverse)
  // Location
  int k_forward;                  // Breakpoint diagonal (forward)
  int k_reverse;                  // Breakpoint diagonal (reverse)
  wf_offset_t offset_forward;     // Offset (forward)
  wf_offset_t offset_reverse;     // Offset (reverse)
  affine2p_matrix_type component; // Component (M/I/D)
} wf_bialign_breakpoint_t;

typedef struct {
  // Wavefronts
  wavefront_aligner_t* alg_forward;    // Forward aligner
  wavefront_aligner_t* alg_reverse;    // Reverse aligner
  wavefront_aligner_t* alg_subsidiary; // Subsidiary aligner
  // Operators
  void (*wf_align_compute)(wavefront_aligner_t* const,const int);
} wavefront_bialigner_t;

/*
 * Setup
 */
wavefront_bialigner_t* wavefront_bialigner_new(
    wavefront_aligner_attr_t* const attributes,
    wavefront_plot_t* const plot);
void wavefront_bialigner_reap(
    wavefront_bialigner_t* const wf_bialigner);
void wavefront_bialigner_delete(
    wavefront_bialigner_t* const wf_bialigner);

/*
 * Sequences
 */
void wavefront_bialigner_set_sequences_ascii(
    wavefront_bialigner_t* const wf_bialigner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length);
void wavefront_bialigner_set_sequences_lambda(
    wavefront_bialigner_t* const wf_bialigner,
    alignment_match_funct_t match_funct,
    void* match_funct_arguments,
    const int pattern_length,
    const int text_length);
void wavefront_bialigner_set_sequences_packed2bits(
    wavefront_bialigner_t* const wf_bialigner,
    const uint8_t* const pattern,
    const int pattern_length,
    const uint8_t* const text,
    const int text_length);
void wavefront_bialigner_set_sequences_bounds(
    wavefront_bialigner_t* const wf_bialigner,
    const int pattern_begin,
    const int pattern_end,
    const int text_begin,
    const int text_end);

/*
 * Accessors
 */
uint64_t wavefront_bialigner_get_size(
    wavefront_bialigner_t* const wf_bialigner);
void wavefront_bialigner_set_heuristic(
    wavefront_bialigner_t* const wf_bialigner,
    wavefront_heuristic_t* const heuristic);
void wavefront_bialigner_set_max_alignment_score(
    wavefront_bialigner_t* const wf_bialigner,
    const int max_alignment_score);
void wavefront_bialigner_set_max_memory(
    wavefront_bialigner_t* const wf_bialigner,
    const uint64_t max_memory_resident,
    const uint64_t max_memory_abort);
void wavefront_bialigner_set_max_num_threads(
    wavefront_bialigner_t* const wf_bialigner,
    const int max_num_threads);
void wavefront_bialigner_set_min_offsets_per_thread(
    wavefront_bialigner_t* const wf_bialigner,
    const int min_offsets_per_thread);
#endif /* WAVEFRONT_BIALIGNER_H_ */
