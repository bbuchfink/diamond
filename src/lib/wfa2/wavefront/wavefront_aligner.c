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

#include "wavefront_aligner.h"
#include "wavefront_components.h"
#include "wavefront_heuristic.h"
#include "wavefront_plot.h"

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
char* wavefront_align_strerror(const int error_code) {
  if (error_code > 0) {
    fprintf(stderr,"[WFA] Internal alignment error code (%d)",error_code);
    exit(1);
  }
  return wf_error_msg[error_code+3];
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
  // Custom matching functions
  wf_aligner->match_funct = attributes->match_funct;
  wf_aligner->match_funct_arguments = attributes->match_funct_arguments;
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
  wf_aligner->sequences = NULL;
  // CIGAR
  const int cigar_length = (score_only) ? 10 : 2*(PATTERN_LENGTH_INIT+TEXT_LENGTH_INIT);
  wf_aligner->cigar = cigar_new(cigar_length,wf_aligner->mm_allocator);
  // System
  wf_aligner->system = attributes->system;
  // Return
  return wf_aligner;
}
void wavefront_aligner_reap(
    wavefront_aligner_t* const wf_aligner) {
  // Padded sequences
  if (wf_aligner->sequences != NULL) {
    strings_padded_delete(wf_aligner->sequences);
    wf_aligner->sequences = NULL;
  }
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
  // Padded sequences
  if (wf_aligner->sequences != NULL) {
    strings_padded_delete(wf_aligner->sequences);
  }
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
    mm_allocator_delete(wf_aligner->mm_allocator);
  }
}
/*
 * Span configuration
 */
void wavefront_aligner_set_alignment_end_to_end(
    wavefront_aligner_t* const wf_aligner) {
  wf_aligner->alignment_form.span = alignment_end2end;
}
void wavefront_aligner_set_alignment_free_ends(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free) {
  wf_aligner->alignment_form.span = alignment_endsfree;
  wf_aligner->alignment_form.pattern_begin_free = pattern_begin_free;
  wf_aligner->alignment_form.pattern_end_free = pattern_end_free;
  wf_aligner->alignment_form.text_begin_free = text_begin_free;
  wf_aligner->alignment_form.text_end_free = text_end_free;
}
/*
 * Heuristic configuration
 */
void wavefront_aligner_set_heuristic_none(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_heuristic_set_none(&wf_aligner->heuristic);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_heuristic_inherit(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
void wavefront_aligner_set_heuristic_banded_static(
    wavefront_aligner_t* const wf_aligner,
    const int band_min_k,
    const int band_max_k) {
  wavefront_heuristic_set_banded_static(&wf_aligner->heuristic,band_min_k,band_max_k);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_heuristic_inherit(wf_aligner->bialigner,&wf_aligner->heuristic);
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
    wavefront_bialigner_heuristic_inherit(wf_aligner->bialigner,&wf_aligner->heuristic);
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
    wavefront_bialigner_heuristic_inherit(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
void wavefront_aligner_set_heuristic_xdrop(
    wavefront_aligner_t* const wf_aligner,
    const int xdrop,
    const int score_steps) {
  wavefront_heuristic_set_xdrop(&wf_aligner->heuristic,xdrop,score_steps);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_heuristic_inherit(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
void wavefront_aligner_set_heuristic_zdrop(
    wavefront_aligner_t* const wf_aligner,
    const int ydrop,
    const int score_steps) {
  wavefront_heuristic_set_zdrop(&wf_aligner->heuristic,ydrop,score_steps);
  if (wf_aligner->bialigner != NULL) {
    wavefront_bialigner_heuristic_inherit(wf_aligner->bialigner,&wf_aligner->heuristic);
  }
}
/*
 * Match-funct configuration
 */
void wavefront_aligner_set_match_funct(
    wavefront_aligner_t* const wf_aligner,
    int (*match_funct)(int,int,void*),
    void* const match_funct_arguments) {
  wf_aligner->match_funct = match_funct;
  wf_aligner->match_funct_arguments = match_funct_arguments;
}
/*
 * System configuration
 */
void wavefront_aligner_set_max_alignment_score(
    wavefront_aligner_t* const wf_aligner,
    const int max_alignment_score) {
  wf_aligner->system.max_alignment_score = max_alignment_score;
}
void wavefront_aligner_set_max_memory(
    wavefront_aligner_t* const wf_aligner,
    const uint64_t max_memory_resident,
    const uint64_t max_memory_abort) {
  wf_aligner->system.max_memory_resident = max_memory_resident;
  wf_aligner->system.max_memory_abort = max_memory_abort;
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
/*
 * Display
 */
void wavefront_aligner_print_type(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner) {
  if (wf_aligner->align_mode_tag == NULL) {
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
  } else {
    fprintf(stream,"%s",wf_aligner->align_mode_tag);
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
void wavefront_aligner_print_mode(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner) {
  fprintf(stream,"(%s,",(wf_aligner->alignment_scope==compute_score)?"Score":"Alg");
  switch (wf_aligner->memory_mode) {
    case wavefront_memory_high: fprintf(stream,"MHigh)"); break;
    case wavefront_memory_med: fprintf(stream,"MMed)"); break;
    case wavefront_memory_low: fprintf(stream,"MLow)"); break;
    case wavefront_memory_ultralow: fprintf(stream,"BiWFA)"); break;
  }
}

