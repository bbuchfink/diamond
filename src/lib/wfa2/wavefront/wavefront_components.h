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
 * DESCRIPTION: WaveFront components module
 */

#ifndef WAVEFRONT_WAVEFRONT_COMPONENTS_H_
#define WAVEFRONT_WAVEFRONT_COMPONENTS_H_

#include "../utils/commons.h"
#include "../wavefront/wavefront.h"
#include "../wavefront/wavefront_backtrace_buffer.h"
#include "../wavefront/wavefront_penalties.h"

/*
 * Wavefront Components
 */
typedef struct {
  // Configuration
  bool memory_modular;                         // Memory strategy (modular wavefronts)
  bool bt_piggyback;                           // Backtrace Piggyback
  // Wavefronts dimensions
  int num_wavefronts;                          // Total number of allocated wavefronts
  int max_score_scope;                         // Maximum score-difference between dependent wavefronts
  int historic_max_hi;                         // Maximum WF hi-limit seen during current alignment
  int historic_min_lo;                         // Minimum WF lo-limit seen during current alignment
  // Wavefronts
  wavefront_t** mwavefronts;                   // M-wavefronts
  wavefront_t** i1wavefronts;                  // I1-wavefronts
  wavefront_t** i2wavefronts;                  // I2-wavefronts
  wavefront_t** d1wavefronts;                  // D1-wavefronts
  wavefront_t** d2wavefronts;                  // D2-wavefronts
  wavefront_t* wavefront_null;                 // Null wavefront (orthogonal reading)
  wavefront_t* wavefront_victim;               // Dummy wavefront (orthogonal writing)
  // BT-Buffer
  wf_backtrace_buffer_t* bt_buffer;            // Backtrace Buffer
  // MM
  mm_allocator_t* mm_allocator;                // MM-Allocator
} wavefront_components_t;

/*
 * Setup
 */
void wavefront_components_allocate(
    wavefront_components_t* const wf_components,
    const int max_pattern_length,
    const int max_text_length,
    wavefront_penalties_t* const penalties,
    const bool memory_modular,
    const bool bt_piggyback,
    mm_allocator_t* const mm_allocator);
void wavefront_components_reap(
    wavefront_components_t* const wf_components);
void wavefront_components_clear(
    wavefront_components_t* const wf_components);
void wavefront_components_free(
    wavefront_components_t* const wf_components);

/*
 * Resize
 */
void wavefront_components_resize(
    wavefront_components_t* const wf_components,
    const int max_pattern_length,
    const int max_text_length,
    wavefront_penalties_t* const penalties);
void wavefront_components_resize_null__victim(
    wavefront_components_t* const wf_components,
    const int lo,
    const int hi);

/*
 * Compact
 */
void wavefront_components_compact_bt_buffer(
    wavefront_components_t* const wf_components,
    const int score,
    const int verbose);

#endif /* WAVEFRONT_WAVEFRONT_COMPONENTS_H_ */
