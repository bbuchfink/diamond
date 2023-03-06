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
 * DESCRIPTION: WaveFront aligner data structure
 */

#include "utils/commons.h"
#include "wavefront_aligner.h"
#include "wavefront_components.h"
#include "wavefront_heuristic.h"
#include "wavefront_plot.h"
#include "wavefront_compute.h"
#include "wavefront_sequences.h"

/*
 * Configuration
 */
#define PATTERN_LENGTH_INIT 1000
#define TEXT_LENGTH_INIT    1000

/*
 * Error messages
 */
char* wf_error_msg[] =
{
  /* WF_STATUS_OOM                  == -3 */ "[WFA] Alignment failed. Maximum memory threshold reached",
  /* WF_STATUS_MAX_SCORE_REACHED    == -2 */ "[WFA] Alignment failed. Maximum score reached",
  /* WF_STATUS_UNFEASIBLE           == -1 */ "[WFA] Alignment unfeasible (possible due to heuristic parameters)",
  /* WF_STATUS_SUCCESSFUL           ==  0 */ "[WFA] Alignment finished successfully",
};
char* wf_error_msg_short[] =
{
  /* WF_STATUS_OOM                  == -3 */ "OOM",
  /* WF_STATUS_MAX_SCORE_REACHED    == -2 */ "MaxScore",
  /* WF_STATUS_UNFEASIBLE           == -1 */ "Unfeasible",
  /* WF_STATUS_SUCCESSFUL           ==  0 */ "OK",
};
char* wavefront_align_strerror(const int error_code) {
  if (error_code > 0) {
    fprintf(stderr,"[WFA] Internal alignment error code (%d)",error_code);
    exit(1);
  }
  return wf_error_msg[error_code+3];
}
char* wavefront_align_strerror_short(const int error_code) {
  return wf_error_msg_short[error_code+3];
}
/*
 * Setup
 */
wavefront_aligner_t* wavefront_aligner_init_mm(
    mm_allocator_t* mm_allocator,
    const bool memory_modular,
    const bool bt_piggyback,
    const bool bi_alignment) {
  // MM
  bool mm_allocator_own;
  if (mm_allocator == NULL) {
    mm_allocator = mm_allocator_new((bi_alignment) ? BUFFER_SIZE_4K : BUFFER_SIZE_4M);
    mm_allocator_own = true;
  } else {
    mm_allocator_own = false;
  }
  // Handler
  wavefront_aligner_t* const wf_aligner =
      mm_allocator_alloc(mm_allocator,wavefront_aligner_t);
  // Configure MM
  wf_aligner->mm_allocator = mm_allocator;
  wf_aligner->mm_allocator_own = mm_allocator_own;
  // Slab
  if (bi_alignment) {
    wf_aligner->wavefront_slab = NULL;
  } else {
    const wf_slab_mode_t slab_mode = (memory_modular) ? wf_slab_reuse : wf_slab_tight;
    wf_aligner->wavefront_slab = wavefront_slab_new(1000,bt_piggyback,slab_mode,wf_aligner->mm_allocator);
  }
  // Return
  return wf_aligner;
}
void wavefront_aligner_init_penalties(
    wavefront_aligner_t* const wf_aligner,
    wavefront_aligner_attr_t* const attributes) {
  switch (attributes->distance_metric) {
    case indel:
      wavefront_penalties_set_indel(&wf_aligner->penalties);
      break;
    case edit:
      wavefront_penalties_set_edit(&wf_aligner->penalties);
      break;
    case gap_linear:
      wavefront_penalties_set_linear(
          &wf_aligner->penalties,
          &attributes->linear_penalties);
      break;
    case gap_affine:
      wavefront_penalties_set_affine(
          &wf_aligner->penalties,
          &attributes->affine_penalties);
      break;
    case gap_affine_2p:
      wavefront_penalties_set_affine2p(
          &wf_aligner->penalties,
          &attributes->affine2p_penalties);
      break;
  }
}
void wavefront_aligner_init_heuristic(
    wavefront_aligner_t* const wf_aligner,
    wavefront_aligner_attr_t* const attributes) {
  // Parameters
  wavefront_heuristic_t* const wf_heuristic = &attributes->heuristic;
  // Select and configure heuristics
  if (wf_heuristic->strategy == wf_heuristic_none) {
    wavefront_heuristic_set_none(&wf_aligner->heuristic);
  } else {
    // Reset
    wf_aligner->heuristic.strategy = 0;
    // WF-Adaptive
    if (wf_heuristic->strategy & wf_heuristic_wfadaptive) {
      wavefront_heuristic_set_wfadaptive(
          &wf_aligner->heuristic,wf_heuristic->min_wavefront_length,
          wf_heuristic->max_distance_threshold,wf_heuristic->steps_between_cutoffs);
    } else if (wf_heuristic->strategy & wf_heuristic_wfmash) {
      wavefront_heuristic_set_wfmash(
          &wf_aligner->heuristic,wf_heuristic->min_wavefront_length,
          wf_heuristic->max_distance_threshold,wf_heuristic->steps_between_cutoffs);
    }
    // Drops
    if (wf_heuristic->strategy & wf_heuristic_xdrop) {
      wavefront_heuristic_set_xdrop(&wf_aligner->heuristic,
          wf_heuristic->xdrop,wf_heuristic->steps_between_cutoffs);
    } else if (wf_heuristic->strategy & wf_heuristic_zdrop) {
      wavefront_heuristic_set_zdrop(&wf_aligner->heuristic,
          wf_heuristic->zdrop,wf_heuristic->steps_between_cutoffs);
    }
    // Banded
    if (wf_heuristic->strategy & wf_heuristic_banded_static) {
      wavefront_heuristic_set_banded_static(&wf_aligner->heuristic,
          wf_heuristic->min_k,wf_heuristic->max_k);
    } else if (wf_heuristic->strategy & wf_heuristic_banded_adaptive) {
      wavefront_heuristic_set_banded_adaptive(&wf_aligner->heuristic,
          wf_heuristic->min_k,wf_heuristic->max_k,wf_heuristic->steps_between_cutoffs);
    }
  }
}
void wavefront_aligner_init_alignment(
    wavefront_aligner_t* const wf_aligner,
    wavefront_aligner_attr_t* const attributes,
    const bool memory_modular,
    const bool bt_piggyback,
    const bool bi_alignment) {
  // Mode
  wf_aligner->align_mode = (bi_alignment) ? wf_align_biwfa : wf_align_regular;
  wf_aligner->align_mode_tag = NULL;
  // Score & form
  wf_aligner->alignment_scope = attributes->alignment_scope;
  wf_aligner->alignment_form = attributes->alignment_form;
  // Penalties
  wavefront_aligner_init_penalties(wf_aligner,attributes);
  // Memory mode
  wf_aligner->memory_mode = attributes->memory_mode;
  wavefront_aligner_init_heuristic(wf_aligner,attributes);
}
wavefront_aligner_t* wavefront_aligner_new(
    wavefront_aligner_attr_t* attributes) {
  // Parameters
  if (attributes == NULL) attributes = &wavefront_aligner_attr_default;
  const bool score_only = (attributes->alignment_scope == compute_score);
  const bool memory_succint =
      attributes->memory_mode == wavefront_memory_med ||
      attributes->memory_mode == wavefront_memory_low;
  const bool memory_modular = score_only || memory_succint;
  const bool bt_piggyback = !score_only && memory_succint;
  const bool bi_alignment = (attributes->memory_mode == wavefront_memory_ultralow);
  // Handler
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_init_mm(
      attributes->mm_allocator,memory_modular,bt_piggyback,bi_alignment);
  // Plot
  if (attributes->plot.enabled) {
    wf_aligner->plot = wavefront_plot_new(attributes->distance_metric,
        PATTERN_LENGTH_INIT,TEXT_LENGTH_INIT,&attributes->plot);
  } else {
    wf_aligner->plot = NULL;
  }
  // Alignment
  wavefront_aligner_init_alignment(wf_aligner,attributes,memory_modular,bt_piggyback,bi_alignment);
  if (bi_alignment) {
    wf_aligner->bialigner = wavefront_bialigner_new(attributes,wf_aligner->plot);
  } else {
    wf_aligner->bialigner = NULL;
    // Wavefront components
    wavefront_components_allocate(
        &wf_aligner->wf_components,PATTERN_LENGTH_INIT,TEXT_LENGTH_INIT,
        &wf_aligner->penalties,memory_modular,bt_piggyback,
        wf_aligner->mm_allocator);
  }
  // Sequences
  wavefront_sequences_allocate(&wf_aligner->sequences);
  // CIGAR
  const int cigar_length = (score_only) ? 10 : 2*(PATTERN_LENGTH_INIT+TEXT_LENGTH_INIT);
  wf_aligner->cigar = cigar_new(cigar_length);
  // System
  wf_aligner->system = attributes->system;
  // Return
  return wf_aligner;
}
void wavefront_aligner_reap(
    wavefront_aligner_t* const wf_aligner) {
  // Select alignment mode
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_reap(wf_aligner->bialigner);
  } else {
    // Wavefront components
    wavefront_components_reap(&wf_aligner->wf_components);
    // Slab
    wavefront_slab_reap(wf_aligner->wavefront_slab);
  }
}
void wavefront_aligner_delete(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  mm_allocator_t* const mm_allocator = wf_aligner->mm_allocator;
  const bool mm_allocator_own = wf_aligner->mm_allocator_own;
  // Sequences
  wavefront_sequences_free(&wf_aligner->sequences);
  // Select alignment mode
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_delete(wf_aligner->bialigner);
  } else {
    // Wavefront components
    wavefront_components_free(&wf_aligner->wf_components);
    // Slab
    wavefront_slab_delete(wf_aligner->wavefront_slab);
  }
  // CIGAR
  cigar_free(wf_aligner->cigar);
  // Plot
  if (wf_aligner->plot != NULL && wf_aligner->align_mode <= 1) {
    wavefront_plot_delete(wf_aligner->plot);
  }
  // MM
  mm_allocator_free(mm_allocator,wf_aligner);
  if (mm_allocator_own) {
    mm_allocator_delete(mm_allocator);
  }
}
/*
 * Initialize wf-alignment conditions
 */
void wavefront_aligner_init_wf_m(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  wavefront_slab_t* const wavefront_slab = wf_aligner->wavefront_slab;
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  alignment_form_t* const form = &wf_aligner->alignment_form;
  // Consider ends-free
  const int hi = (penalties->match==0) ? form->text_begin_free : 0;
  const int lo = (penalties->match==0) ? -form->pattern_begin_free : 0;
  // Compute dimensions
  int effective_lo, effective_hi;
  wavefront_compute_limits_output(wf_aligner,lo,hi,&effective_lo,&effective_hi);
  // Initialize end2end (wavefront zero)
  wf_components->mwavefronts[0] = wavefront_slab_allocate(wavefront_slab,effective_lo,effective_hi);
  wf_components->mwavefronts[0]->offsets[0] = 0;
  wf_components->mwavefronts[0]->lo = lo;
  wf_components->mwavefronts[0]->hi = hi;
  // Store initial BT-piggypack element
  if (wf_components->bt_piggyback) {
    const bt_block_idx_t block_idx = wf_backtrace_buffer_init_block(wf_components->bt_buffer,0,0);
    wf_components->mwavefronts[0]->bt_pcigar[0] = 0;
    wf_components->mwavefronts[0]->bt_prev[0] = block_idx;
  }
  // Initialize ends-free
  if (form->span == alignment_endsfree && penalties->match == 0) {
    // Text begin-free
    const int text_begin_free = form->text_begin_free;
    int h;
    for (h=1;h<=text_begin_free;++h) {
      const int k = DPMATRIX_DIAGONAL(h,0);
      wf_components->mwavefronts[0]->offsets[k] = DPMATRIX_OFFSET(h,0);
      if (wf_components->bt_piggyback) {
        const bt_block_idx_t block_idx = wf_backtrace_buffer_init_block(wf_components->bt_buffer,0,h);
        wf_components->mwavefronts[0]->bt_pcigar[k] = 0;
        wf_components->mwavefronts[0]->bt_prev[k] = block_idx;
      }
    }
    // Pattern begin-free
    const int pattern_begin_free = form->pattern_begin_free;
    int v;
    for (v=1;v<=pattern_begin_free;++v) {
      const int k = DPMATRIX_DIAGONAL(0,v);
      wf_components->mwavefronts[0]->offsets[k] = DPMATRIX_OFFSET(0,v);
      if (wf_components->bt_piggyback) {
        const bt_block_idx_t block_idx = wf_backtrace_buffer_init_block(wf_components->bt_buffer,v,0);
        wf_components->mwavefronts[0]->bt_pcigar[k] = 0;
        wf_components->mwavefronts[0]->bt_prev[k] = block_idx;
      }
    }
  }
  // Nullify unused WFs
  if (distance_metric <= gap_linear) return;
  wf_components->d1wavefronts[0] = NULL;
  wf_components->i1wavefronts[0] = NULL;
  if (distance_metric==gap_affine) return;
  wf_components->d2wavefronts[0] = NULL;
  wf_components->i2wavefronts[0] = NULL;
}
void wavefront_aligner_init_wf(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  wavefront_slab_t* const wavefront_slab = wf_aligner->wavefront_slab;
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Init wavefronts
  if (wf_aligner->component_begin == affine2p_matrix_M) {
    // Initialize
    wavefront_aligner_init_wf_m(wf_aligner);
    // Nullify unused WFs
    if (distance_metric <= gap_linear) return;
    wf_components->i1wavefronts[0] = NULL;
    wf_components->d1wavefronts[0] = NULL;
    if (distance_metric == gap_affine) return;
    wf_components->i2wavefronts[0] = NULL;
    wf_components->d2wavefronts[0] = NULL;
  } else {
    // Compute dimensions
    int effective_lo, effective_hi; // Effective lo/hi
    wavefront_compute_limits_output(wf_aligner,0,0,&effective_lo,&effective_hi);
    wavefront_t* const wavefront = wavefront_slab_allocate(wavefront_slab,effective_lo,effective_hi);
    // Initialize
    switch (wf_aligner->component_begin) {
      case affine2p_matrix_I1:
        wf_components->mwavefronts[0] = NULL;
        wf_components->i1wavefronts[0] = wavefront;
        wf_components->i1wavefronts[0]->offsets[0] = 0;
        wf_components->i1wavefronts[0]->lo = 0;
        wf_components->i1wavefronts[0]->hi = 0;
        wf_components->d1wavefronts[0] = NULL;
        // Nullify unused WFs
        if (distance_metric==gap_affine) return;
        wf_components->i2wavefronts[0] = NULL;
        wf_components->d2wavefronts[0] = NULL;
        break;
      case affine2p_matrix_I2:
        wf_components->mwavefronts[0] = NULL;
        wf_components->i1wavefronts[0] = NULL;
        wf_components->d1wavefronts[0] = NULL;
        wf_components->i2wavefronts[0] = wavefront;
        wf_components->i2wavefronts[0]->offsets[0] = 0;
        wf_components->i2wavefronts[0]->lo = 0;
        wf_components->i2wavefronts[0]->hi = 0;
        wf_components->d2wavefronts[0] = NULL;
        break;
      case affine2p_matrix_D1:
        wf_components->mwavefronts[0] = NULL;
        wf_components->i1wavefronts[0] = NULL;
        wf_components->d1wavefronts[0] = wavefront;
        wf_components->d1wavefronts[0]->offsets[0] = 0;
        wf_components->d1wavefronts[0]->lo = 0;
        wf_components->d1wavefronts[0]->hi = 0;
        // Nullify unused WFs
        if (distance_metric==gap_affine) return;
        wf_components->i2wavefronts[0] = NULL;
        wf_components->d2wavefronts[0] = NULL;
        break;
      case affine2p_matrix_D2:
        wf_components->mwavefronts[0] = NULL;
        wf_components->i1wavefronts[0] = NULL;
        wf_components->d1wavefronts[0] = NULL;
        wf_components->i2wavefronts[0] = NULL;
        wf_components->d2wavefronts[0] = wavefront;
        wf_components->d2wavefronts[0]->offsets[0] = 0;
        wf_components->d2wavefronts[0]->lo = 0;
        wf_components->d2wavefronts[0]->hi = 0;
        break;
      default:
        break;
    }
  }
}
/*
 * Span configuration
 */
void wavefront_aligner_set_alignment_end_to_end(
    wavefront_aligner_t* const wf_aligner) {
  wf_aligner->alignment_form.span = alignment_end2end;
  wf_aligner->alignment_form.extension = false;
}
void wavefront_aligner_set_alignment_free_ends(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free) {
  wf_aligner->alignment_form.span = alignment_endsfree;
  wf_aligner->alignment_form.extension = false;
  wf_aligner->alignment_form.pattern_begin_free = pattern_begin_free;
  wf_aligner->alignment_form.pattern_end_free = pattern_end_free;
  wf_aligner->alignment_form.text_begin_free = text_begin_free;
  wf_aligner->alignment_form.text_end_free = text_end_free;
}
void wavefront_aligner_set_alignment_extension(
    wavefront_aligner_t* const wf_aligner) {
  wf_aligner->alignment_form.span = alignment_endsfree;
  wf_aligner->alignment_form.extension = true;
}
/*
 * Heuristic configuration
 */
void wavefront_aligner_set_heuristic_none(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_heuristic_set_none(&wf_aligner->heuristic);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_heuristic(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
void wavefront_aligner_set_heuristic_banded_static(
    wavefront_aligner_t* const wf_aligner,
    const int band_min_k,
    const int band_max_k) {
  wavefront_heuristic_set_banded_static(&wf_aligner->heuristic,band_min_k,band_max_k);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_heuristic(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
void wavefront_aligner_set_heuristic_banded_adaptive(
    wavefront_aligner_t* const wf_aligner,
    const int band_min_k,
    const int band_max_k,
    const int score_steps) {
  wavefront_heuristic_set_banded_adaptive(
      &wf_aligner->heuristic,band_min_k,band_max_k,score_steps);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_heuristic(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
void wavefront_aligner_set_heuristic_wfadaptive(
    wavefront_aligner_t* const wf_aligner,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int score_steps) {
  wavefront_heuristic_set_wfadaptive(
      &wf_aligner->heuristic,
      min_wavefront_length,max_distance_threshold,score_steps);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_heuristic(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
void wavefront_aligner_set_heuristic_wfmash(
    wavefront_aligner_t* const wf_aligner,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int score_steps) {
  wavefront_heuristic_set_wfmash(
      &wf_aligner->heuristic,
      min_wavefront_length,max_distance_threshold,score_steps);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_heuristic(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
void wavefront_aligner_set_heuristic_xdrop(
    wavefront_aligner_t* const wf_aligner,
    const int xdrop,
    const int score_steps) {
  wavefront_heuristic_set_xdrop(&wf_aligner->heuristic,xdrop,score_steps);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_heuristic(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
void wavefront_aligner_set_heuristic_zdrop(
    wavefront_aligner_t* const wf_aligner,
    const int ydrop,
    const int score_steps) {
  wavefront_heuristic_set_zdrop(&wf_aligner->heuristic,ydrop,score_steps);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_heuristic(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
/*
 * System configuration
 */
void wavefront_aligner_set_max_alignment_score(
    wavefront_aligner_t* const wf_aligner,
    const int max_alignment_score) {
  wf_aligner->system.max_alignment_score = max_alignment_score;
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_max_alignment_score(
        wf_aligner->bialigner,max_alignment_score);
  }
}
void wavefront_aligner_set_max_memory(
    wavefront_aligner_t* const wf_aligner,
    const uint64_t max_memory_resident,
    const uint64_t max_memory_abort) {
  wf_aligner->system.max_memory_resident = max_memory_resident;
  wf_aligner->system.max_memory_abort = max_memory_abort;
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_max_memory(
        wf_aligner->bialigner,max_memory_resident,max_memory_abort);
  }
}
void wavefront_aligner_set_max_num_threads(
    wavefront_aligner_t* const wf_aligner,
    const int max_num_threads) {
  wf_aligner->system.max_num_threads = max_num_threads;
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_max_num_threads(
        wf_aligner->bialigner,max_num_threads);
  }
}
void wavefront_aligner_set_min_offsets_per_thread(
    wavefront_aligner_t* const wf_aligner,
    const int min_offsets_per_thread) {
  wf_aligner->system.min_offsets_per_thread = min_offsets_per_thread;
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_set_min_offsets_per_thread(
        wf_aligner->bialigner,min_offsets_per_thread);
  }
}
/*
 * Utils
 */
uint64_t wavefront_aligner_get_size(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  // Bialigner
  uint64_t sub_aligners = 0;
  if (wf_aligner->bialigner != NULL) {
    return wavefront_bialigner_get_size(wf_aligner->bialigner);
  } else {
    // Compute aligner size
    const uint64_t bt_buffer_size = (wf_components->bt_buffer) ?
        wf_backtrace_buffer_get_size_allocated(wf_components->bt_buffer) : 0;
    const uint64_t slab_size = wavefront_slab_get_size(wf_aligner->wavefront_slab);
    // Return overall size
    return sub_aligners + bt_buffer_size + slab_size;
  }
}
void wavefront_aligner_maxtrim_cigar(
    wavefront_aligner_t* const wf_aligner) {
  switch (wf_aligner->penalties.distance_metric) {
    case indel:
    case edit:
      // Maxtrim does not apply to edit/indel distances
      break;
    case gap_linear:
      cigar_maxtrim_gap_linear(wf_aligner->cigar,&wf_aligner->penalties.linear_penalties);
      break;
    case gap_affine:
      cigar_maxtrim_gap_affine(wf_aligner->cigar,&wf_aligner->penalties.affine_penalties);
      break;
    case gap_affine_2p:
      cigar_maxtrim_gap_affine2p(wf_aligner->cigar,&wf_aligner->penalties.affine2p_penalties);
      break;
  }
}
/*
 * Display
 */
void wavefront_aligner_print_mode(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner) {
  if (wf_aligner->align_mode_tag != NULL) {
    fprintf(stream,"%s::",wf_aligner->align_mode_tag);
  }
  switch (wf_aligner->align_mode) {
    case wf_align_biwfa:
      fprintf(stream,"BiWFA");
      break;
    case wf_align_biwfa_breakpoint_forward:
      fprintf(stream,"BiWFA::Forward");
      break;
    case wf_align_biwfa_breakpoint_reverse:
      fprintf(stream,"BiWFA::Reverse");
      break;
    case wf_align_biwfa_subsidiary:
      fprintf(stream,"BiWFA::SubWFA");
      break;
    default:
      fprintf(stream,"WFA");
      break;
  }
}
void wavefront_aligner_print_scope(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner) {
  const char* const scope_label =
      (wf_aligner->alignment_scope == compute_score) ? "score" : "alignment";
  if (wf_aligner->alignment_form.span == alignment_end2end) {
    fprintf(stream,"(%s,end2end)",scope_label);
  } else {
    fprintf(stream,"(%s,endsfree,%d,%d,%d,%d)",
        scope_label,
        wf_aligner->alignment_form.pattern_begin_free,
        wf_aligner->alignment_form.pattern_end_free,
        wf_aligner->alignment_form.text_begin_free,
        wf_aligner->alignment_form.text_end_free);
  }
}
void wavefront_aligner_print_conf(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner) {
  fprintf(stream,"(");
  switch (wf_aligner->memory_mode) {
    case wavefront_memory_high: fprintf(stream,"MHigh"); break;
    case wavefront_memory_med: fprintf(stream,"MMed"); break;
    case wavefront_memory_low: fprintf(stream,"MLow"); break;
    case wavefront_memory_ultralow: fprintf(stream,"BiWFA"); break;
  }
  if (wf_aligner->system.max_alignment_score == INT_MAX) {
    fprintf(stream,",inf)");
  } else {
    fprintf(stream,",%d)",wf_aligner->system.max_alignment_score);
  }
}
