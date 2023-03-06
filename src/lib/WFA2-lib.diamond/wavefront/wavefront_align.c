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
 * DESCRIPTION: WaveFront alignment module for sequence pairwise alignment
 */

#include "utils/commons.h"
#include "wavefront_align.h"
#include "wavefront_unialign.h"
#include "wavefront_bialign.h"
#include "wavefront_compute.h"
#include "wavefront_compute_edit.h"
#include "wavefront_compute_linear.h"
#include "wavefront_compute_affine.h"
#include "wavefront_compute_affine2p.h"
#include "wavefront_extend.h"
#include "wavefront_backtrace.h"
#include "wavefront_debug.h"

/*
 * Checks
 */
void wavefront_align_checks(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_length,
    const int text_length) {
  alignment_form_t* const form = &wf_aligner->alignment_form;
  if (wf_aligner->bialigner != NULL) {
    const bool ends_free =
        form->pattern_begin_free > 0 ||
        form->pattern_end_free > 0 ||
        form->text_begin_free > 0 ||
        form->text_end_free > 0;
    if (ends_free) {
      fprintf(stderr,"[WFA] BiWFA ends-free has not been tested properly yet (let me know and I'll do it)\n");
      exit(1);
    }
    if (wf_aligner->alignment_form.extension) {
      fprintf(stderr,"[WFA] BiWFA extension is not implemented yet (let me know and I'll add it)\n");
      exit(1);
    }
  }
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  const bool is_heuristic_drop =
      (wf_aligner->heuristic.strategy & wf_heuristic_xdrop) ||
      (wf_aligner->heuristic.strategy & wf_heuristic_zdrop);
  if (is_heuristic_drop && (distance_metric==edit || distance_metric==indel)) {
    fprintf(stderr,"[WFA] Heuristics drops are not compatible with 'edit'/'indel' distance metrics\n");
    exit(1);
  }
  if (form->span == alignment_endsfree) {
    if (form->pattern_begin_free > pattern_length ||
        form->pattern_end_free > pattern_length ||
        form->text_begin_free > text_length ||
        form->text_end_free > text_length) {
      fprintf(stderr,"[WFA] Ends-free parameters must be not larger than the sequences "
          "(P0=%d,Pf=%d,T0=%d,Tf=%d). Must be (P0<=|P|,Pf<=|P|,T0<=|T|,Tf<=|T|) where (|P|,|T|)=(%d,%d)\n",
          form->pattern_begin_free,form->pattern_end_free,
          form->text_begin_free,form->text_end_free,
          pattern_length,text_length);
      exit(1);
    }
  }
}
/*
 * Wavefront Alignment Unidirectional
 */
void wavefront_align_unidirectional_cleanup(
    wavefront_aligner_t* const wf_aligner) {
  // Compute memory used
  uint64_t memory_used = wavefront_aligner_get_size(wf_aligner);
  wf_aligner->align_status.memory_used = memory_used;
  // Reap memory (controlled reaping)
  if (memory_used > wf_aligner->system.max_memory_resident) {
    // Wavefront components
    wavefront_components_reap(&wf_aligner->wf_components);
    // Check memory
    memory_used = wavefront_aligner_get_size(wf_aligner);
    wf_aligner->align_status.memory_used = memory_used;
    // Slab
    if (memory_used > wf_aligner->system.max_memory_resident) {
      wavefront_slab_reap(wf_aligner->wavefront_slab);
      if (wf_aligner->bialigner != NULL) {
        wavefront_bialigner_reap(wf_aligner->bialigner);
      }
    }
  }
}
void wavefront_align_unidirectional(
    wavefront_aligner_t* const wf_aligner) {
  // Wavefront align sequences
  wavefront_unialign_init(wf_aligner,affine2p_matrix_M,affine2p_matrix_M); // Init
  wavefront_unialign(wf_aligner); // Align
  // Finish
  if (wf_aligner->align_status.status == WF_STATUS_MAX_SCORE_REACHED) return; // Alignment paused
  wavefront_align_unidirectional_cleanup(wf_aligner);
}
/*
 * Wavefront Alignment Bidirectional
 */
void wavefront_align_bidirectional(
    wavefront_aligner_t* const wf_aligner) {
  // Bidirectional alignment
  wavefront_bialign(wf_aligner); // Align
  // Finish
  wf_aligner->align_status.memory_used = wavefront_aligner_get_size(wf_aligner);
}
/*
 * Wavefront Alignment Dispatcher
 */
int wavefront_align_lambda(
    wavefront_aligner_t* const wf_aligner,
    alignment_match_funct_t match_funct,
    void* match_funct_arguments,
    const int pattern_length,
    const int text_length) {
  // Checks
  wavefront_align_checks(wf_aligner,pattern_length,text_length);
  wavefront_debug_begin(wf_aligner);
  // Plot
  if (wf_aligner->plot != NULL) wavefront_plot_resize(wf_aligner->plot,pattern_length,text_length);
  // Dispatcher
  if (wf_aligner->bialigner == NULL) {
    // Prepare Sequences
    wavefront_sequences_init_lambda(&wf_aligner->sequences,
        match_funct,match_funct_arguments,
        pattern_length,text_length,false);
    wavefront_align_unidirectional(wf_aligner);
  } else {
    // Prepare Sequences
    wavefront_bialigner_set_sequences_lambda(wf_aligner->bialigner,
        match_funct,match_funct_arguments,
        pattern_length,text_length);
    // Align
    wavefront_align_bidirectional(wf_aligner);
  }
  // DEBUG
  wavefront_debug_end(wf_aligner);
  wavefront_debug_check_correct(wf_aligner);
  // Return
  return wf_aligner->align_status.status;
}
int wavefront_align_packed2bits(
    wavefront_aligner_t* const wf_aligner,
    const uint8_t* const pattern,
    const int pattern_length,
    const uint8_t* const text,
    const int text_length) {
  // Checks
  wavefront_align_checks(wf_aligner,pattern_length,text_length);
  wavefront_debug_begin(wf_aligner);
  // Plot
  if (wf_aligner->plot != NULL) wavefront_plot_resize(wf_aligner->plot,pattern_length,text_length);
  // Dispatcher
  if (wf_aligner->bialigner == NULL) {
    // Prepare Sequences
    wavefront_sequences_init_packed2bits(&wf_aligner->sequences,
        pattern,pattern_length,text,text_length,false);
    wavefront_align_unidirectional(wf_aligner);
  } else {
    // Prepare Sequences
    wavefront_bialigner_set_sequences_packed2bits(wf_aligner->bialigner,
        pattern,pattern_length,text,text_length);
    // Align
    wavefront_align_bidirectional(wf_aligner);
  }
  // DEBUG
  wavefront_debug_end(wf_aligner);
  wavefront_debug_check_correct(wf_aligner);
  // Return
  return wf_aligner->align_status.status;
}
int wavefront_align(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  // Checks
  wavefront_align_checks(wf_aligner,pattern_length,text_length);
  wavefront_debug_begin(wf_aligner);
  // Plot
  if (wf_aligner->plot != NULL) wavefront_plot_resize(wf_aligner->plot,pattern_length,text_length);
  // Dispatcher
  if (wf_aligner->bialigner == NULL) {
    // Prepare Sequences
    wavefront_sequences_init_ascii(&wf_aligner->sequences,
        pattern,pattern_length,text,text_length,false);
    wavefront_align_unidirectional(wf_aligner);
  } else {
    // Prepare Sequences
    wavefront_bialigner_set_sequences_ascii(wf_aligner->bialigner,
        pattern,pattern_length,text,text_length);
    // Align
    wavefront_align_bidirectional(wf_aligner);
  }
  // DEBUG
  wavefront_debug_end(wf_aligner);
  wavefront_debug_check_correct(wf_aligner);
  // Return
  return wf_aligner->align_status.status;
}
/*
 * Wavefront Alignment Resume (Experimental)
 */
int wavefront_align_resume(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  wavefront_align_status_t* const align_status = &wf_aligner->align_status;
  // Check current alignment status
  if (align_status->status != WF_STATUS_MAX_SCORE_REACHED ||
      wf_aligner->bialigner != NULL) {
    fprintf(stderr,"[WFA] Alignment cannot be resumed\n");
    exit(1);
  }
  // Resume aligning sequences
  wavefront_unialign(wf_aligner);
  // Finish alignment
  if (align_status->status == WF_STATUS_MAX_SCORE_REACHED) {
    return WF_STATUS_MAX_SCORE_REACHED; // Alignment paused
  }
  wavefront_align_unidirectional_cleanup(wf_aligner);
  // DEBUG
  wavefront_debug_check_correct(wf_aligner);
  // Return
  return align_status->status;
}
