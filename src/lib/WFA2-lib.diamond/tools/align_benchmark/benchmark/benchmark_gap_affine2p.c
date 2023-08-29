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
 * DESCRIPTION: Gap-affine 2-pieces alignment algorithms wrapper (including WFA)
 */

#include "benchmark/benchmark_gap_affine2p.h"
#include "benchmark/benchmark_check.h"
#include "gap_affine2p/affine2p_matrix.h"
#include "gap_affine2p/affine2p_dp.h"
#include "wavefront/wavefront_align.h"

/*
 * Benchmark Gap-Affine 2-Pieces
 */
void benchmark_gap_affine2p_dp(
    align_input_t* const align_input,
    affine2p_penalties_t* const penalties) {
  // Allocate
  affine2p_matrix_t matrix;
  affine2p_matrix_allocate(
      &matrix,align_input->pattern_length+1,
      align_input->text_length+1,align_input->mm_allocator);
  cigar_t* const cigar = cigar_new(
      align_input->pattern_length+align_input->text_length);
  // Align
  timer_start(&align_input->timer);
  affine2p_dp_align(&matrix,penalties,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,cigar);
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,cigar);
  }
  // Output
  if (align_input->output_file) {
    benchmark_print_output(align_input,gap_affine_2p,false,cigar);
  }
  // Free
  affine2p_matrix_free(&matrix,align_input->mm_allocator);
  cigar_free(cigar);
}
void benchmark_gap_affine2p_wavefront(
    align_input_t* const align_input,
    affine2p_penalties_t* const penalties) {
  // Parameters
  wavefront_aligner_t* const wf_aligner = align_input->wf_aligner;
  // Align
  timer_start(&align_input->timer);
  if (align_input->wfa_match_funct == NULL) {
    wavefront_align(wf_aligner,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length);
  } else {
    wavefront_align_lambda(wf_aligner,
        align_input->wfa_match_funct,align_input->wfa_match_funct_arguments,
        align_input->pattern_length,align_input->text_length);
  }
  timer_stop(&align_input->timer);
  // DEBUG
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,wf_aligner->cigar);
  }
  // Output
  if (align_input->output_file) {
    const int score_only = (wf_aligner->alignment_scope == compute_score);
    benchmark_print_output(align_input,gap_affine_2p,score_only,wf_aligner->cigar);
  }
}
