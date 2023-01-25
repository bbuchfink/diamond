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

#ifndef WAVEFRONT_SLAB_H_
#define WAVEFRONT_SLAB_H_

#include "../utils/commons.h"
#include "../utils/vector.h"
#include "../system/mm_allocator.h"
#include "wavefront.h"

/*
 * Memory Manager for Wavefront
 */
typedef enum {
  wf_slab_reuse = 1, // Keep all wavefronts (Reap only by demand)
  wf_slab_tight = 2, // Reap all if wavefronts are resized
} wf_slab_mode_t;
typedef struct {
  // Attributes
  bool allocate_backtrace;         // WFs require BT-vector
  wf_slab_mode_t slab_mode;        // Slab strategy
  // Wavefront Slabs
  int init_wf_length;              // Initial wf-elements allocated
  int current_wf_length;           // Current wf-elements allocated
  vector_t* wavefronts;            // All wavefronts (wavefront_t*)
  vector_t* wavefronts_free;       // Free wavefronts (wavefront_t*)
  // Stats
  uint64_t memory_used;            // Memory used (Bytes)
  // MM
  mm_allocator_t* mm_allocator;    // MM-Allocator
} wavefront_slab_t;

/*
 * Setup
 */
wavefront_slab_t* wavefront_slab_new(
    const int init_wf_length,
    const bool allocate_backtrace,
    const wf_slab_mode_t slab_mode,
    mm_allocator_t* const mm_allocator);
void wavefront_slab_reap(
    wavefront_slab_t* const wavefront_slab);
void wavefront_slab_clear(
    wavefront_slab_t* const wavefront_slab);
void wavefront_slab_delete(
    wavefront_slab_t* const wavefront_slab);

/*
 * Accessors
 */
void wavefront_slab_set_mode(
    wavefront_slab_t* const wavefront_slab,
    const wf_slab_mode_t slab_mode);

/*
 * Allocator
 */
wavefront_t* wavefront_slab_allocate(
    wavefront_slab_t* const wavefront_slab,
    const int min_lo,
    const int max_hi);
void wavefront_slab_free(
    wavefront_slab_t* const wavefront_slab,
    wavefront_t* const wavefront);

/*
 * Utils
 */
uint64_t wavefront_slab_get_size(
    wavefront_slab_t* const wavefront_slab);

#endif /* WAVEFRONT_SLAB_H_ */


