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
 * DESCRIPTION: Individual WaveFront data structure
 */

#ifndef WAVEFRONT_H_
#define WAVEFRONT_H_

#include "../utils/commons.h"
#include "../system/mm_allocator.h"
#include "wavefront_offset.h"
#include "wavefront_backtrace_buffer.h"

/*
 * Alignment position
 */
typedef struct {
  int score;          // Score
  int k;              // Diagonal
  wf_offset_t offset; // Offset
} wavefront_pos_t;

/*
 * Wavefront
 */
typedef enum {
  wavefront_status_free,
  wavefront_status_busy,
  wavefront_status_deallocated,
} wavefront_status_type;
typedef struct {
  // Dimensions
  bool null;                           // Is null interval?
  int lo;                              // Lowest diagonal (inclusive)
  int hi;                              // Highest diagonal (inclusive)
  // Wavefront elements
  wf_offset_t* offsets;                // Offsets (k-centered)
  wf_offset_t* offsets_mem;            // Offsets base memory (Internal)
  // Piggyback backtrace
  int bt_occupancy_max;                // Maximum number of pcigar-ops stored on the Backtrace-block
  pcigar_t* bt_pcigar;                 // Backtrace-block pcigar (k-centered)
  bt_block_idx_t* bt_prev;             // Backtrace-block previous-index (k-centered)
  pcigar_t* bt_pcigar_mem;             // Backtrace-block (base memory - Internal)
  bt_block_idx_t* bt_prev_mem;         // Backtrace-block previous-index (base memory - Internal)
  // Slab internals
  wavefront_status_type status;        // Wavefront status (memory state)
  int wf_elements_allocated;           // Total wf-elements allocated (max. wf. size)
  int wf_elements_allocated_min;       // Minimum diagonal-element wf-element allocated
  int wf_elements_allocated_max;       // Maximum diagonal-element wf-element allocated
  int wf_elements_init_min;            // Minimum diagonal-element initialized (inclusive)
  int wf_elements_init_max;            // Maximum diagonal-element initialized (inclusive)
} wavefront_t;

/*
 * Wavefront Set
 */
typedef struct {
  /* In Wavefronts*/
  wavefront_t* in_mwavefront_misms;
  wavefront_t* in_mwavefront_open1;
  wavefront_t* in_mwavefront_open2;
  wavefront_t* in_i1wavefront_ext;
  wavefront_t* in_i2wavefront_ext;
  wavefront_t* in_d1wavefront_ext;
  wavefront_t* in_d2wavefront_ext;
  /* Out Wavefronts */
  wavefront_t* out_mwavefront;
  wavefront_t* out_i1wavefront;
  wavefront_t* out_i2wavefront;
  wavefront_t* out_d1wavefront;
  wavefront_t* out_d2wavefront;
} wavefront_set_t;

/*
 * Setup
 */
void wavefront_allocate(
    wavefront_t* const wavefront,
    const int wf_elements_allocated,
    const bool allocate_backtrace,
    mm_allocator_t* const mm_allocator);
void wavefront_resize(
    wavefront_t* const wavefront,
    const int wf_elements_allocated,
    mm_allocator_t* const mm_allocator);
void wavefront_free(
    wavefront_t* const wavefront,
    mm_allocator_t* const mm_allocator);

/*
 * Initialization
 */
void wavefront_init(
    wavefront_t* const wavefront,
    const int min_lo,
    const int max_hi);
void wavefront_init_null(
    wavefront_t* const wavefront,
    const int min_lo,
    const int max_hi);
void wavefront_init_victim(
    wavefront_t* const wavefront,
    const int min_lo,
    const int max_hi);

/*
 * Accessors
 */
void wavefront_set_limits(
    wavefront_t* const wavefront,
    const int lo,
    const int hi);

/*
 * Utils
 */
uint64_t wavefront_get_size(
    wavefront_t* const wavefront);

#endif /* WAVEFRONT_H_ */
