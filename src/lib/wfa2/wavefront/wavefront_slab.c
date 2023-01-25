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
 * DESCRIPTION: WaveFront Slab for fast pre-allocated wavefronts' memory handling
 */

#include "wavefront_slab.h"

/*
 * Constants
 */
#define WF_SLAB_EXPAND_FACTOR          1.5f
#define WF_SLAB_QUEUES_LENGTH_INIT     100

/*
 * Setup
 */
wavefront_slab_t* wavefront_slab_new(
    const int init_wf_length,
    const bool allocate_backtrace,
    const wf_slab_mode_t slab_mode,
    mm_allocator_t* const mm_allocator) {
  // Allocate
  wavefront_slab_t* const wavefront_slab =
      mm_allocator_alloc(mm_allocator,wavefront_slab_t);
  // Attributes
  wavefront_slab->allocate_backtrace = allocate_backtrace;
  wavefront_slab->slab_mode = slab_mode;
  // Wavefront Slabs
  wavefront_slab->init_wf_length = init_wf_length;
  wavefront_slab->current_wf_length = init_wf_length;
  wavefront_slab->wavefronts = vector_new(WF_SLAB_QUEUES_LENGTH_INIT,wavefront_t*);
  wavefront_slab->wavefronts_free = vector_new(WF_SLAB_QUEUES_LENGTH_INIT,wavefront_t*);
  // Stats
  wavefront_slab->memory_used = 0;
  // MM
  wavefront_slab->mm_allocator = mm_allocator;
  // Return
  return wavefront_slab;
}
void wavefront_slab_reap_free(
    wavefront_slab_t* const wavefront_slab) {
  // Parameters
  wavefront_t** const wavefronts = vector_get_mem(wavefront_slab->wavefronts,wavefront_t*);
  const int num_wavefronts = vector_get_used(wavefront_slab->wavefronts);
  mm_allocator_t* const mm_allocator = wavefront_slab->mm_allocator;
  // Remove "deallocated" and free wavefronts
  int i, valid_idx = 0;
  for (i=0;i<num_wavefronts;++i) {
    switch (wavefronts[i]->status) {
      case wavefront_status_deallocated:
        mm_allocator_free(mm_allocator,wavefronts[i]); // Delete handler
        break;
      case wavefront_status_busy:
        wavefronts[valid_idx++] = wavefronts[i]; // Valid wavefront
        break;
      case wavefront_status_free:
        wavefront_free(wavefronts[i],mm_allocator); // Free wavefront
        wavefront_slab->memory_used -= wavefront_get_size(wavefronts[i]);
        mm_allocator_free(mm_allocator,wavefronts[i]); // Delete handler
        break;
    }
  }
  vector_set_used(wavefront_slab->wavefronts,valid_idx);
  vector_clear(wavefront_slab->wavefronts_free);
}
void wavefront_slab_reap_repurpose(
    wavefront_slab_t* const wavefront_slab) {
  // Parameters
  const int current_wf_length = wavefront_slab->current_wf_length;
  wavefront_t** const wavefronts = vector_get_mem(wavefront_slab->wavefronts,wavefront_t*);
  const int num_wavefronts = vector_get_used(wavefront_slab->wavefronts);
  mm_allocator_t* const mm_allocator = wavefront_slab->mm_allocator;
  // Clear free
  vector_reserve(wavefront_slab->wavefronts_free,num_wavefronts,false);
  wavefront_t** const wavefronts_free = vector_get_mem(wavefront_slab->wavefronts_free,wavefront_t*);
  // Remove "deallocated" and repurpose all we can of current wf-length
  int i, valid_idx = 0;
  for (i=0;i<num_wavefronts;++i) {
    switch (wavefronts[i]->status) {
      case wavefront_status_deallocated:
        mm_allocator_free(mm_allocator,wavefronts[i]); // Delete handler
        break;
      case wavefront_status_busy:
      case wavefront_status_free:
        if (wavefronts[i]->wf_elements_allocated == current_wf_length) {
          wavefronts[i]->status = wavefront_status_free; // Set free
          wavefronts[valid_idx] = wavefronts[i]; // Valid wavefront
          wavefronts_free[valid_idx] = wavefronts[i]; // Free wavefront
          valid_idx++;
        } else {
          wavefront_free(wavefronts[i],mm_allocator); // Free wavefront
          wavefront_slab->memory_used -= wavefront_get_size(wavefronts[i]);
          mm_allocator_free(mm_allocator,wavefronts[i]); // Delete handler
        }
        break;
    }
  }
  vector_set_used(wavefront_slab->wavefronts,valid_idx);
  vector_set_used(wavefront_slab->wavefronts_free,valid_idx);
}
void wavefront_slab_reap(
    wavefront_slab_t* const wavefront_slab) {
  // Back to initial size
  wavefront_slab->current_wf_length = wavefront_slab->init_wf_length;
  wavefront_slab_reap_repurpose(wavefront_slab); // Repurpose all wavefronts
}
void wavefront_slab_clear(
    wavefront_slab_t* const wavefront_slab) {
  // Select slab mode
  switch (wavefront_slab->slab_mode) {
    case wf_slab_reuse:
      wavefront_slab_reap_repurpose(wavefront_slab);
      break;
    case wf_slab_tight:
      // Back to initial size
      wavefront_slab->current_wf_length = wavefront_slab->init_wf_length;
      wavefront_slab_reap_repurpose(wavefront_slab);
      break;
  }
}
void wavefront_slab_delete(
    wavefront_slab_t* const wavefront_slab) {
  // Parameters
  mm_allocator_t* const mm_allocator = wavefront_slab->mm_allocator;
  // Delete free vector
  vector_delete(wavefront_slab->wavefronts_free);
  // Free wavefronts
  wavefront_t** const wavefronts =
      vector_get_mem(wavefront_slab->wavefronts,wavefront_t*);
  const int num_wavefronts = vector_get_used(wavefront_slab->wavefronts);
  int i;
  for (i=0;i<num_wavefronts;++i) {
    if (wavefronts[i]->status == wavefront_status_deallocated) {
      mm_allocator_free(mm_allocator,wavefronts[i]); // Delete handler
    } else {
      wavefront_free(wavefronts[i],mm_allocator); // Free wavefront
      mm_allocator_free(mm_allocator,wavefronts[i]); // Delete handler
    }
  }
  vector_delete(wavefront_slab->wavefronts);
  // Handler
  mm_allocator_free(wavefront_slab->mm_allocator,wavefront_slab);
}
/*
 * Accessors
 */
void wavefront_slab_set_mode(
    wavefront_slab_t* const wavefront_slab,
    const wf_slab_mode_t slab_mode) {
  // Check mode
  if (slab_mode != wavefront_slab->slab_mode) {
    // Change mode
    wavefront_slab->slab_mode = slab_mode;
    // Reap
    wavefront_slab->current_wf_length = wavefront_slab->init_wf_length;
    wavefront_slab_reap_repurpose(wavefront_slab);
  }
}
/*
 * Slab Allocator
 */
wavefront_t* wavefront_slab_allocate_new(
    wavefront_slab_t* const wavefront_slab,
    const int wf_length_requested,
    const int min_lo,
    const int max_hi) {
  // Allocate a new wavefront
  mm_allocator_t* const mm_allocator = wavefront_slab->mm_allocator;
  wavefront_t* const wavefront = mm_allocator_alloc(mm_allocator,wavefront_t);
  wavefront_allocate(wavefront,wf_length_requested,wavefront_slab->allocate_backtrace,mm_allocator);
  vector_insert(wavefront_slab->wavefronts,wavefront,wavefront_t*);
  wavefront_slab->memory_used += wavefront_get_size(wavefront);
  // Init wavefront
  wavefront->status = wavefront_status_busy;
  wavefront_init(wavefront,min_lo,max_hi);
  // Return
  return wavefront;
}
wavefront_t* wavefront_slab_allocate_free(
    wavefront_slab_t* const wavefront_slab,
    const int min_lo,
    const int max_hi) {
  // Parameters
  vector_t* const wavefronts_free = wavefront_slab->wavefronts_free;
  // Reuse wavefront
  wavefront_t* const wavefront = *(vector_get_last_elm(wavefronts_free,wavefront_t*));
  vector_dec_used(wavefronts_free);
  // Init wavefront
  wavefront->status = wavefront_status_busy;
  wavefront_init(wavefront,min_lo,max_hi);
  // Return
  return wavefront;
}
wavefront_t* wavefront_slab_allocate(
    wavefront_slab_t* const wavefront_slab,
    const int min_lo,
    const int max_hi) {
  // Parameters
  vector_t* const wavefronts_free = wavefront_slab->wavefronts_free;
  const int wf_length = WAVEFRONT_LENGTH(min_lo,max_hi);
  // Check slab-mode
  if (wavefront_slab->slab_mode == wf_slab_reuse) {
    // Check max-length of pre-allocated wavefronts
    if (wf_length > wavefront_slab->current_wf_length) {
      const int proposed_wf_length = (float)wf_length * WF_SLAB_EXPAND_FACTOR;
      wavefront_slab->current_wf_length = proposed_wf_length; // New slab size
      wavefront_slab_reap_free(wavefront_slab); // Reap free wavefronts
    }
    // Check for a free wavefront (pre-allocated in the slab)
    if (vector_get_used(wavefronts_free) > 0) {
      return wavefront_slab_allocate_free(wavefront_slab,min_lo,max_hi);
    } else {
      // Allocate a new wavefront
      return wavefront_slab_allocate_new(wavefront_slab,
          wavefront_slab->current_wf_length,min_lo,max_hi);
    }
  } else { // wf_slab_tight
    if (wf_length <= wavefront_slab->init_wf_length) {
      // Check for a free wavefront (pre-allocated in the slab)
      if (vector_get_used(wavefronts_free) > 0) {
        return wavefront_slab_allocate_free(wavefront_slab,min_lo,max_hi);
      } else {
        return wavefront_slab_allocate_new(wavefront_slab,
            wavefront_slab->init_wf_length,min_lo,max_hi); // Allocate new
      }
    } else {
      return wavefront_slab_allocate_new(wavefront_slab,
          wf_length,min_lo,max_hi); // Allocate new
    }
  }
}
void wavefront_slab_free(
    wavefront_slab_t* const wavefront_slab,
    wavefront_t* const wavefront) {
  // Check reasons to repurpose wavefront (NOTE: Tight-mode never slab_frees())
  //   (A) Reuse-mode and wavefront has current wf-length
  //   (B) Tight-mode and wavefront has init wf-length
  const int wf_length = wavefront->wf_elements_allocated;
  const bool repurpose_reuse =
      (wavefront_slab->slab_mode == wf_slab_reuse) &&
      (wf_length == wavefront_slab->current_wf_length);
  const bool repurpose_tight =
      (wavefront_slab->slab_mode == wf_slab_tight) &&
      (wf_length == wavefront_slab->init_wf_length);
  if (repurpose_reuse || repurpose_tight) {
    // Return wavefront to slab as free (Good job recycling)
    wavefront->status = wavefront_status_free;
    vector_insert(wavefront_slab->wavefronts_free,wavefront,wavefront_t*);
  } else {
    // Delete wavefront
    wavefront->status = wavefront_status_deallocated;
    wavefront_slab->memory_used -= wavefront_get_size(wavefront);
    wavefront_free(wavefront,wavefront_slab->mm_allocator);
  }
}
/*
 * Utils
 */
uint64_t wavefront_slab_get_size(
    wavefront_slab_t* const wavefront_slab) {
  return wavefront_slab->memory_used;
}



