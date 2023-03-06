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
 * DESCRIPTION: Benchmark utils
 */

#include "benchmark_check.h"
#include "alignment/score_matrix.h"
#include "indel/indel_dp.h"
#include "edit/edit_dp.h"
#include "gap_linear/nw.h"
#include "gap_affine/affine_matrix.h"
#include "gap_affine/swg.h"
#include "gap_affine2p/affine2p_dp.h"

/*
 * Checker (given the correct CIGAR+SCORE)
 */
void benchmark_check_alignment_using_solution(
    align_input_t* const align_input,
    cigar_t* const cigar_computed,
    const int score_computed,
    cigar_t* const cigar_correct,
    const int score_correct) {
  counter_add(&(align_input->align),1);
  counter_add(&(align_input->align_score_total),ABS(score_computed));
  // Debug
  if (align_input->debug_flags) {
    // Display info
    if (align_input->debug_flags & ALIGN_DEBUG_DISPLAY_INFO) {
      benchmark_print_alignment(stderr,align_input,score_computed,cigar_computed,-1,NULL);
    }
    // Check correct
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_CORRECT) {
      bool correct = cigar_check_alignment(stderr,
          align_input->pattern,align_input->pattern_length,
          align_input->text,align_input->text_length,
          cigar_computed,align_input->verbose);
      if (!correct) {
        // Print
        if (align_input->verbose) {
          fprintf(stderr,"INCORRECT ALIGNMENT\n");
          benchmark_print_alignment(stderr,align_input,-1,cigar_computed,-1,NULL);
        }
        // Quit
        return;
      } else {
        counter_add(&(align_input->align_correct),1);
      }
      // CIGAR Stats
      int i;
      counter_add(&(align_input->align_bases),align_input->pattern_length);
      for (i=cigar_computed->begin_offset;i<cigar_computed->end_offset;++i) {
        switch (cigar_computed->operations[i]) {
          case 'M': counter_add(&(align_input->align_matches),1); break;
          case 'X': counter_add(&(align_input->align_mismatches),1); break;
          case 'I': counter_add(&(align_input->align_ins),1); break;
          case 'D': default: counter_add(&(align_input->align_del),1); break;
        }
      }
    }
    // Check score
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_SCORE) {
      if (score_computed != score_correct) {
        // Print
        if (align_input->verbose) {
          benchmark_print_alignment(
              stderr,align_input,
              score_computed,cigar_computed,
              score_correct,cigar_correct);
          fprintf(stderr,"(#%d)\t INACCURATE SCORE computed=%d\tcorrect=%d\n",
              align_input->sequence_id,score_computed,score_correct);
        }
        counter_add(&(align_input->align_score_diff),ABS(score_computed-score_correct));
        // Quit
        return;
      } else {
        counter_add(&(align_input->align_score),1);
      }
    }
    // Check alignment
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_ALIGNMENT) {
      if (cigar_cmp(cigar_computed,cigar_correct) != 0) {
        // Print
        if (align_input->verbose) {
          fprintf(stderr,"INACCURATE ALIGNMENT\n");
          benchmark_print_alignment(
              stderr,align_input,
              -1,cigar_computed,
              -1,cigar_correct);
        }
        // Quit
        return;
      } else {
        counter_add(&(align_input->align_cigar),1);
      }
    }
  }
}
/*
 * Specialized Checkers
 */
void benchmark_check_alignment_indel(
    align_input_t* const align_input,
    cigar_t* const cigar_computed) {
  score_matrix_t score_matrix;
  score_matrix_allocate(
      &score_matrix,align_input->pattern_length+1,
      align_input->text_length+1,align_input->mm_allocator);
  cigar_t* const cigar = cigar_new(
      align_input->pattern_length+align_input->text_length);
  indel_dp_compute(&score_matrix,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,cigar);
  const int score_correct = cigar_score_edit(cigar);
  const int score_computed = cigar_score_edit(cigar_computed);
  // Check alignment
  benchmark_check_alignment_using_solution(
      align_input,cigar_computed,score_computed,
      cigar,score_correct);
  // Free
  score_matrix_free(&score_matrix);
  cigar_free(cigar);
}
void benchmark_check_alignment_edit(
    align_input_t* const align_input,
    cigar_t* const cigar_computed) {
  score_matrix_t score_matrix;
  score_matrix_allocate(
      &score_matrix,align_input->pattern_length+1,
      align_input->text_length+1,align_input->mm_allocator);
  cigar_t* const cigar = cigar_new(
      align_input->pattern_length+align_input->text_length);
  if (align_input->check_bandwidth <= 0) {
    edit_dp_align(&score_matrix,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length,cigar);
  } else {
    edit_dp_align_banded(&score_matrix,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length,
        align_input->check_bandwidth,cigar);
  }
  const int score_correct = cigar_score_edit(cigar);
  const int score_computed = cigar_score_edit(cigar_computed);
  // Check alignment
  benchmark_check_alignment_using_solution(
      align_input,cigar_computed,score_computed,
      cigar,score_correct);
  // Free
  score_matrix_free(&score_matrix);
  cigar_free(cigar);
}
void benchmark_check_alignment_gap_linear(
    align_input_t* const align_input,
    cigar_t* const cigar_computed) {
  // Compute correct
  score_matrix_t score_matrix;
  score_matrix_allocate(
      &score_matrix,align_input->pattern_length+1,
      align_input->text_length+1,align_input->mm_allocator);
  cigar_t* const cigar = cigar_new(
      align_input->pattern_length+align_input->text_length);
  nw_align(&score_matrix,
      align_input->check_linear_penalties,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,cigar);
  const int score_correct = cigar_score_gap_linear(
      cigar,align_input->check_linear_penalties);
  const int score_computed = cigar_score_gap_linear(
      cigar_computed,align_input->check_linear_penalties);
  // Check alignment
  benchmark_check_alignment_using_solution(
      align_input,cigar_computed,score_computed,
      cigar,score_correct);
  // Free
  score_matrix_free(&score_matrix);
  cigar_free(cigar);
}
void benchmark_check_alignment_gap_affine(
    align_input_t* const align_input,
    cigar_t* const cigar_computed) {
  // Allocate
  affine_matrix_t affine_matrix;
  affine_matrix_allocate(
      &affine_matrix,align_input->pattern_length+1,
      align_input->text_length+1,align_input->mm_allocator);
  cigar_t* const cigar = cigar_new(
      align_input->pattern_length+align_input->text_length);
  // Compute correct
  if (align_input->check_bandwidth <= 0) {
    swg_align(&affine_matrix,align_input->check_affine_penalties,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length,cigar);
  } else {
    swg_align_banded(&affine_matrix,align_input->check_affine_penalties,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length,
        align_input->check_bandwidth,cigar);
  }
  const int score_correct = cigar_score_gap_affine(
      cigar,align_input->check_affine_penalties);
  const int score_computed = cigar_score_gap_affine(
      cigar_computed,align_input->check_affine_penalties);
  // Check alignment
  benchmark_check_alignment_using_solution(
      align_input,cigar_computed,score_computed,
      cigar,score_correct);
  // Free
  affine_matrix_free(&affine_matrix,align_input->mm_allocator);
  cigar_free(cigar);
}
void benchmark_check_alignment_gap_affine2p(
    align_input_t* const align_input,
    cigar_t* const cigar_computed) {
  // Allocate
  affine2p_matrix_t affine_matrix;
  affine2p_matrix_allocate(
      &affine_matrix,align_input->pattern_length+1,
      align_input->text_length+1,align_input->mm_allocator);
  cigar_t* const cigar = cigar_new(
      align_input->pattern_length+align_input->text_length);
  // Compute correct
  affine2p_dp_align(
      &affine_matrix,align_input->check_affine2p_penalties,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,cigar);
  const int score_correct = cigar_score_gap_affine2p(
      cigar,align_input->check_affine2p_penalties);
  const int score_computed = cigar_score_gap_affine2p(
      cigar_computed,align_input->check_affine2p_penalties);
  // Check alignment
  benchmark_check_alignment_using_solution(
      align_input,cigar_computed,score_computed,
      cigar,score_correct);
  // Free
  affine2p_matrix_free(&affine_matrix,align_input->mm_allocator);
  cigar_free(cigar);
}
/*
 * Check
 */
void benchmark_check_alignment(
    align_input_t* const align_input,
    cigar_t* const cigar_computed) {
  // Compute correct CIGAR
  if ((align_input->debug_flags & ALIGN_DEBUG_CHECK_SCORE) ||
      (align_input->debug_flags & ALIGN_DEBUG_CHECK_ALIGNMENT)) {
    // Check gap-affine-2p
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE2P) {
      benchmark_check_alignment_gap_affine2p(align_input,cigar_computed);
    }
    // Check gap-affine
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE) {
      benchmark_check_alignment_gap_affine(align_input,cigar_computed);
    }
    // Check gap-linear
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_LINEAR) {
      benchmark_check_alignment_gap_linear(align_input,cigar_computed);
    }
    // Check edit
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_DISTANCE_METRIC_EDIT) {
      benchmark_check_alignment_edit(align_input,cigar_computed);
    }
    // Check indel
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_DISTANCE_METRIC_INDEL) {
      benchmark_check_alignment_indel(align_input,cigar_computed);
    }
  } else if ((align_input->debug_flags & ALIGN_DEBUG_CHECK_CORRECT) ||
             (align_input->debug_flags & ALIGN_DEBUG_DISPLAY_INFO)) { // Check CORRECT or just DISPLAY
    // Delegate check alignment
    benchmark_check_alignment_using_solution(
        align_input,cigar_computed,-1,NULL,-1);
  }
}

