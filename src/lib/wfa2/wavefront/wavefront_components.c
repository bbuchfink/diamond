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
 * DESCRIPTION: WaveFront aligner components
 */

#include "wavefront_components.h"
#include "../utils/bitmap.h"
#include "../system/profiler_timer.h"

/*
 * Configuration
 */
#define WF_NULL_INIT_LO     (-1024)
#define WF_NULL_INIT_HI     ( 1024)
#define WF_NULL_INIT_LENGTH WAVEFRONT_LENGTH(WF_NULL_INIT_LO,WF_NULL_INIT_HI)

/*
 * Compute dimensions
 */
void wavefront_components_dimensions_edit(
    wavefront_components_t* const wf_components,
    const int max_pattern_length,
    const int max_text_length,
    int* const max_score_scope,
    int* const num_wavefronts) {
  // Compute max-scope
  *max_score_scope = 2;
  // Dimensions
  if (wf_components->memory_modular) {
    *num_wavefronts = 2;
  } else {
    *num_wavefronts = MAX(max_pattern_length,max_text_length);
  }
}
void wavefront_components_dimensions_linear(
    wavefront_components_t* const wf_components,
    wavefront_penalties_t* const penalties,
    const int max_pattern_length,
    const int max_text_length,
    int* const max_score_scope,
    int* const num_wavefronts) {
  // Compute max-scope
  *max_score_scope = MAX(penalties->mismatch,penalties->gap_opening1) + 1;
  // Dimensions
  if (wf_components->memory_modular) {
    *num_wavefronts = *max_score_scope;
  } else {
    const int abs_seq_diff = ABS(max_pattern_length-max_text_length);
    const int max_score_misms = MIN(max_pattern_length,max_text_length) * penalties->mismatch;
    const int max_score_indel = penalties->gap_opening1 * abs_seq_diff;
    *num_wavefronts = max_score_misms + max_score_indel;
  }
}
void wavefront_components_dimensions_affine(
    wavefront_components_t* const wf_components,
    wavefront_penalties_t* const penalties,
    const int max_pattern_length,
    const int max_text_length,
    int* const max_score_scope,
    int* const num_wavefronts) {
  // Compute max-scope
  const int max_score_scope_indel = penalties->gap_opening1+penalties->gap_extension1;
  *max_score_scope = MAX(max_score_scope_indel,penalties->mismatch) + 1;
  // Dimensions
  if (wf_components->memory_modular) {
    *num_wavefronts = *max_score_scope;
  } else {
    const int abs_seq_diff = ABS(max_pattern_length-max_text_length);
    const int max_score_misms = MIN(max_pattern_length,max_text_length) * penalties->mismatch;
    const int max_score_indel = penalties->gap_opening1 + abs_seq_diff * penalties->gap_extension1;
    *num_wavefronts = max_score_misms + max_score_indel;
  }
}
void wavefront_components_dimensions_affine2p(
    wavefront_components_t* const wf_components,
    wavefront_penalties_t* const penalties,
    const int max_pattern_length,
    const int max_text_length,
    int* const max_score_scope,
    int* const num_wavefronts) {
  // Compute max-scope
  const int max_score_scope_indel =
      MAX(penalties->gap_opening1+penalties->gap_extension1,
          penalties->gap_opening2+penalties->gap_extension2);
  *max_score_scope = MAX(max_score_scope_indel,penalties->mismatch) + 1;
  // Dimensions
  if (wf_components->memory_modular) {
    *num_wavefronts = *max_score_scope;
  } else {
    const int abs_seq_diff = ABS(max_pattern_length-max_text_length);
    const int max_score_misms = MIN(max_pattern_length,max_text_length) * penalties->mismatch;
    const int max_score_indel1 = penalties->gap_opening1 + abs_seq_diff * penalties->gap_extension1;
    const int max_score_indel2 = penalties->gap_opening2 + abs_seq_diff * penalties->gap_extension2;
    const int max_score_indel = MIN(max_score_indel1,max_score_indel2);
    *num_wavefronts = max_score_misms + max_score_indel;
  }
}
void wavefront_components_dimensions(
    wavefront_components_t* const wf_components,
    wavefront_penalties_t* const penalties,
    const int max_pattern_length,
    const int max_text_length,
    int* const max_score_scope,
    int* const num_wavefronts) {
  // Switch attending to distance-metric
  switch (penalties->distance_metric) {
    case indel:
    case edit:
      wavefront_components_dimensions_edit(
          wf_components,
          max_pattern_length,max_text_length,
          max_score_scope,num_wavefronts);
      break;
    case gap_linear:
      wavefront_components_dimensions_linear(
          wf_components,penalties,
          max_pattern_length,max_text_length,
          max_score_scope,num_wavefronts);
      break;
    case gap_affine:
      wavefront_components_dimensions_affine(
          wf_components,penalties,
          max_pattern_length,max_text_length,
          max_score_scope,num_wavefronts);
      break;
    case gap_affine_2p:
      wavefront_components_dimensions_affine2p(
          wf_components,penalties,
          max_pattern_length,max_text_length,
          max_score_scope,num_wavefronts);
      break;
  }
  // Clear historic
  wf_components->historic_max_hi = 0;
  wf_components->historic_min_lo = 0;
}
/*
 * Setup
 */
void wavefront_components_allocate_wf(
    wavefront_components_t* const wf_components,
    const int max_pattern_length,
    const int max_text_length,
    const distance_metric_t distance_metric) {
  // Parameters
  const int num_wavefronts = wf_components->num_wavefronts;
  const bool init_wf = wf_components->memory_modular;
  mm_allocator_t* const mm_allocator = wf_components->mm_allocator;
  // Allocate wavefronts
  wf_components->mwavefronts = mm_allocator_calloc(mm_allocator,num_wavefronts,wavefront_t*,init_wf);
  if (distance_metric <= gap_linear) {
    wf_components->i1wavefronts = NULL;
    wf_components->d1wavefronts = NULL;
    wf_components->i2wavefronts = NULL;
    wf_components->d2wavefronts = NULL;
  } else {
    wf_components->i1wavefronts = mm_allocator_calloc(mm_allocator,num_wavefronts,wavefront_t*,init_wf);
    wf_components->d1wavefronts = mm_allocator_calloc(mm_allocator,num_wavefronts,wavefront_t*,init_wf);
    if (distance_metric == gap_affine) {
      wf_components->i2wavefronts = NULL;
      wf_components->d2wavefronts = NULL;
    } else {
      wf_components->i2wavefronts = mm_allocator_calloc(mm_allocator,num_wavefronts,wavefront_t*,init_wf);
      wf_components->d2wavefronts = mm_allocator_calloc(mm_allocator,num_wavefronts,wavefront_t*,init_wf);
    }
  }
}
void wavefront_components_allocate(
    wavefront_components_t* const wf_components,
    const int max_pattern_length,
    const int max_text_length,
    wavefront_penalties_t* const penalties,
    const bool memory_modular,
    const bool bt_piggyback,
    mm_allocator_t* const mm_allocator) {
  // Configuration
  wf_components->memory_modular = memory_modular;
  wf_components->bt_piggyback = bt_piggyback;
  wf_components->mm_allocator = mm_allocator; // MM
  // Allocate wavefronts
  wavefront_components_dimensions(
      wf_components,penalties,
      max_pattern_length,max_text_length,
      &wf_components->max_score_scope,
      &wf_components->num_wavefronts);
  wavefront_components_allocate_wf(wf_components,
      max_pattern_length,max_text_length,penalties->distance_metric);
  // Allocate victim wavefront (outside slab)
  wavefront_t* const wavefront_victim = mm_allocator_alloc(mm_allocator,wavefront_t);
  wavefront_allocate(wavefront_victim,WF_NULL_INIT_LENGTH,bt_piggyback,mm_allocator);
  wavefront_init_victim(wavefront_victim,WF_NULL_INIT_LO,WF_NULL_INIT_HI);
  wf_components->wavefront_victim = wavefront_victim;
  // Allocate null wavefront (outside slab)
  wavefront_t* const wavefront_null = mm_allocator_alloc(mm_allocator,wavefront_t);
  wavefront_allocate(wavefront_null,WF_NULL_INIT_LENGTH,bt_piggyback,mm_allocator);
  wavefront_init_null(wavefront_null,WF_NULL_INIT_LO,WF_NULL_INIT_HI);
  wf_components->wavefront_null = wavefront_null;
  // BT-Buffer
  wf_components->bt_buffer = (bt_piggyback) ? wf_backtrace_buffer_new(mm_allocator) : NULL;
}
void wavefront_components_reap(
    wavefront_components_t* const wf_components) {
  // BT-Buffer
  if (wf_components->bt_buffer) wf_backtrace_buffer_reap(wf_components->bt_buffer);
}
void wavefront_components_clear(
    wavefront_components_t* const wf_components) {
  // Wavefronts components
  if (wf_components->memory_modular) {
    const int num_wavefronts = wf_components->num_wavefronts;
    const int wf_size = num_wavefronts*sizeof(wavefront_t*);
    memset(wf_components->mwavefronts,0,wf_size);
    if (wf_components->i1wavefronts) memset(wf_components->i1wavefronts,0,wf_size);
    if (wf_components->d1wavefronts) memset(wf_components->d1wavefronts,0,wf_size);
    if (wf_components->i2wavefronts) memset(wf_components->i2wavefronts,0,wf_size);
    if (wf_components->d2wavefronts) memset(wf_components->d2wavefronts,0,wf_size);
  }
  wf_components->historic_max_hi = 0;
  wf_components->historic_min_lo = 0;
  // BT-Buffer
  if (wf_components->bt_buffer) wf_backtrace_buffer_clear(wf_components->bt_buffer);
}
void wavefront_components_free_wf(
    wavefront_components_t* const wf_components) {
  // Parameters
  mm_allocator_t* const mm_allocator = wf_components->mm_allocator;
  // Wavefronts components
  mm_allocator_free(mm_allocator,wf_components->mwavefronts);
  if (wf_components->i1wavefronts) mm_allocator_free(mm_allocator,wf_components->i1wavefronts);
  if (wf_components->d1wavefronts) mm_allocator_free(mm_allocator,wf_components->d1wavefronts);
  if (wf_components->i2wavefronts) mm_allocator_free(mm_allocator,wf_components->i2wavefronts);
  if (wf_components->d2wavefronts) mm_allocator_free(mm_allocator,wf_components->d2wavefronts);
}
void wavefront_components_free(
    wavefront_components_t* const wf_components) {
  // Parameters
  mm_allocator_t* const mm_allocator = wf_components->mm_allocator;
  // Wavefronts components
  wavefront_components_free_wf(wf_components);
  // Null wavefront
  wavefront_free(wf_components->wavefront_null,mm_allocator);
  mm_allocator_free(mm_allocator,wf_components->wavefront_null);
  // Victim wavefront
  wavefront_free(wf_components->wavefront_victim,mm_allocator);
  mm_allocator_free(mm_allocator,wf_components->wavefront_victim);
  // BT-Buffer
  if (wf_components->bt_buffer) wf_backtrace_buffer_delete(wf_components->bt_buffer);
}
/*
 * Resize
 */
void wavefront_components_resize(
    wavefront_components_t* const wf_components,
    const int max_pattern_length,
    const int max_text_length,
    wavefront_penalties_t* const penalties) {
  // Compute dimensions
  int num_wavefronts = 0;
  wavefront_components_dimensions(
      wf_components,penalties,
      max_pattern_length,max_text_length,
      &wf_components->max_score_scope,&num_wavefronts);
  // Resize wavefronts components (if needed)
  if (num_wavefronts > wf_components->num_wavefronts) {
    wf_components->num_wavefronts = num_wavefronts;
    wavefront_components_free_wf(wf_components);
    wavefront_components_allocate_wf(wf_components,
        max_pattern_length,max_text_length,penalties->distance_metric);
    // BT-Buffer
    if (wf_components->bt_buffer) wf_backtrace_buffer_clear(wf_components->bt_buffer);
  } else {
    wavefront_components_clear(wf_components);
  }
}
void wavefront_components_resize_null__victim(
    wavefront_components_t* const wf_components,
    const int lo,
    const int hi) {
  // Resize null/victim wavefronts (if needed)
  if (lo-1 < wf_components->wavefront_null->wf_elements_init_min ||
      hi+1 > wf_components->wavefront_null->wf_elements_init_max) {
    // Parameters
    mm_allocator_t* const mm_allocator = wf_components->mm_allocator;
    // Expand and leave some leeway
    const int wf_inc = (WAVEFRONT_LENGTH(lo,hi)*3)/2;
    const int proposed_lo = lo - wf_inc/2;
    const int proposed_hi = hi + wf_inc/2;
    const int proposed_wavefront_length = WAVEFRONT_LENGTH(proposed_lo,proposed_hi);
    // Reallocate victim wavefront
    wavefront_resize(wf_components->wavefront_victim,proposed_wavefront_length,mm_allocator);
    wavefront_init_victim(wf_components->wavefront_victim,proposed_lo,proposed_hi);
    // Allocate null wavefront
    wavefront_resize(wf_components->wavefront_null,proposed_wavefront_length,mm_allocator);
    wavefront_init_null(wf_components->wavefront_null,proposed_lo,proposed_hi);
  }
}
/*
 * Mark wavefronts
 */
void wavefront_components_mark_backtrace(
    wf_backtrace_buffer_t* const bt_buffer,
    bitmap_t* const bitmap,
    wavefront_t* const wavefront) {
  // Parameters
  wf_offset_t* const offsets = wavefront->offsets;
  bt_block_idx_t* const bt_prev = wavefront->bt_prev;
  const int lo = wavefront->lo;
  const int hi = wavefront->hi;
  // Mark all wavefront backtraces (batch mode)
  wf_backtrace_buffer_mark_backtrace_batch(
      bt_buffer,offsets+lo,bt_prev+lo,hi-lo+1,bitmap);
  //  int i;
  //  for (i=lo;i<=hi;++i) {
  //    if (offsets[i]>=0) wf_backtrace_buffer_mark_backtrace(bt_buffer,bt_prev[i],bitmap);
  //  }
}
void wavefront_components_mark_wavefronts(
    wavefront_components_t* const wf_components,
    bitmap_t* const bitmap,
    const int score) {
  // Parameters
  wf_backtrace_buffer_t* const bt_buffer = wf_components->bt_buffer;
  const int max_score_scope = wf_components->max_score_scope;
  // Mark Active Working Set (AWS)
  int i;
  for (i=0;i<max_score_scope;++i) {
    // Compute score
    const int score_mod = (score-i) % wf_components->max_score_scope;
    // Mark M-wavefront
    wavefront_t* const mwavefront = wf_components->mwavefronts[score_mod];
    if (mwavefront!=NULL) wavefront_components_mark_backtrace(bt_buffer,bitmap,mwavefront);
    // Mark (I1/D1)-wavefronts
    if (wf_components->i1wavefronts != NULL) {
      wavefront_t* const i1wavefront = wf_components->i1wavefronts[score_mod];
      if (i1wavefront!=NULL) wavefront_components_mark_backtrace(bt_buffer,bitmap,i1wavefront);
      wavefront_t* const d1wavefront = wf_components->d1wavefronts[score_mod];
      if (d1wavefront!=NULL) wavefront_components_mark_backtrace(bt_buffer,bitmap,d1wavefront);
      // Mark (I2/D2)-wavefronts
      if (wf_components->i2wavefronts != NULL) {
        wavefront_t* const i2wavefront = wf_components->i2wavefronts[score_mod];
        if (i2wavefront!=NULL) wavefront_components_mark_backtrace(bt_buffer,bitmap,i2wavefront);
        wavefront_t* const d2wavefront = wf_components->d2wavefronts[score_mod];
        if (d2wavefront!=NULL) wavefront_components_mark_backtrace(bt_buffer,bitmap,d2wavefront);
      }
    }
  }
  // Update counters in marked bitmap
  bitmap_update_counters(bitmap);
}
/*
 * Translate block-idxs
 */
void wavefront_components_translate_idx(
    wavefront_components_t* const wf_components,
    bitmap_t* const bitmap,
    wavefront_t* const wavefront) {
  // Parameters
  wf_offset_t* const offsets = wavefront->offsets;
  bt_block_idx_t* const bt_prev = wavefront->bt_prev;
  const int lo = wavefront->lo;
  const int hi = wavefront->hi;
  const bt_block_idx_t num_compacted_blocks = wf_components->bt_buffer->num_compacted_blocks;
  // Translate all wavefront block-idxs
  int k;
  for (k=lo;k<=hi;++k) {
    if (offsets[k]>=0) {  // NOTE bt_prev[k] >= num_compacted_blocks
      bt_prev[k] = (bt_prev[k]==BT_BLOCK_IDX_NULL) ?
          BT_BLOCK_IDX_NULL :
          num_compacted_blocks + bitmap_erank(bitmap,bt_prev[k]);
    }
  }
}
void wavefront_components_translate_wavefronts(
    wavefront_components_t* const wf_components,
    bitmap_t* const bitmap,
    const int score) {
  // Mark Active Working Set (AWS)
  const int max_score_scope = wf_components->max_score_scope;
  int i;
  for (i=0;i<max_score_scope;++i) {
    // Compute score
    const int score_mod = (score-i) % wf_components->max_score_scope;
    // Mark M-wavefront
    wavefront_t* const mwavefront = wf_components->mwavefronts[score_mod];
    if (mwavefront!=NULL) wavefront_components_translate_idx(wf_components,bitmap,mwavefront);
    // Mark (I1/D1)-wavefronts
    if (wf_components->i1wavefronts != NULL) {
      wavefront_t* const i1wavefront = wf_components->i1wavefronts[score_mod];
      if (i1wavefront!=NULL) wavefront_components_translate_idx(wf_components,bitmap,i1wavefront);
      wavefront_t* const d1wavefront = wf_components->d1wavefronts[score_mod];
      if (d1wavefront!=NULL) wavefront_components_translate_idx(wf_components,bitmap,d1wavefront);
      // Mark (I2/D2)-wavefronts
      if (wf_components->i2wavefronts != NULL) {
        wavefront_t* const i2wavefront = wf_components->i2wavefronts[score_mod];
        if (i2wavefront!=NULL) wavefront_components_translate_idx(wf_components,bitmap,i2wavefront);
        wavefront_t* const d2wavefront = wf_components->d2wavefronts[score_mod];
        if (d2wavefront!=NULL) wavefront_components_translate_idx(wf_components,bitmap,d2wavefront);
      }
    }
  }
}
/*
 * Compact
 */
void wavefront_components_compact_bt_buffer(
    wavefront_components_t* const wf_components,
    const int score,
    const int verbose) {
  // PROFILE
  profiler_timer_t timer;
  if (verbose >= 3) { timer_reset(&timer); timer_start(&timer); }
  // Parameters
  wf_backtrace_buffer_t* const bt_buffer = wf_components->bt_buffer;
  const uint64_t bt_buffer_used = wf_backtrace_buffer_get_used(bt_buffer);
  // Allocate bitmap
  bitmap_t* const bitmap = bitmap_new(bt_buffer_used,wf_components->mm_allocator);
  // Mark Active Working Set (AWS)
  wavefront_components_mark_wavefronts(wf_components,bitmap,score);
  // Compact marked blocks (also translates idxs to compacted positions)
  const bt_block_idx_t total_compacted_blocks =
      wf_backtrace_buffer_compact_marked(bt_buffer,bitmap,verbose);
  // Translate Active Working Set (AWS)
  wavefront_components_translate_wavefronts(wf_components,bitmap,score);
  // Set new compacted blocks
  wf_backtrace_buffer_set_num_compacted_blocks(bt_buffer,total_compacted_blocks);
  // Free
  bitmap_delete(bitmap);
  // PROFILE
  if (verbose >= 3) {
    timer_stop(&timer);
    fprintf(stderr,"[");
    timer_print_total(stderr,&timer);
    fprintf(stderr,"]\n");
  }
}

