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
 * DESCRIPTION: WaveFront aligner data structure attributes
 */

#ifndef WAVEFRONT_ATTRIBUTES_H_
#define WAVEFRONT_ATTRIBUTES_H_

#include "../utils/commons.h"
#include "../alignment/cigar.h"
#include "../alignment/affine_penalties.h"
#include "../alignment/affine2p_penalties.h"
#include "../alignment/linear_penalties.h"
#include "../system/profiler_timer.h"
#include "../system/mm_allocator.h"

#include "wavefront_penalties.h"
#include "wavefront_plot.h"
#include "wavefront_display.h"
#include "wavefront_heuristic.h"

/*
 * Alignment scope
 */
typedef enum {
  compute_score,           // Only distance/score
  compute_alignment,       // Full alignment CIGAR
} alignment_scope_t;
typedef enum {
  alignment_end2end,       // End-to-end alignment (aka global)
  alignment_endsfree,      // Ends-free alignment  (semiglobal, glocal, etc)
} alignment_span_t;
typedef struct {
  // Mode
  alignment_span_t span;   // Alignment form (End-to-end/Ends-free)
  // Ends-free
  int pattern_begin_free;  // Allow free-gap at the beginning of the pattern
  int pattern_end_free;    // Allow free-gap at the end of the pattern
  int text_begin_free;     // Allow free-gap at the beginning of the text
  int text_end_free;       // Allow free-gap at the end of the text
} alignment_form_t;

/*
 * Custom extend-match function, e.g.:
 *
 *   typedef struct {
 *     char* pattern;
 *     int pattern_length;
 *     char* text;
 *     int text_length;
 *   } match_function_params_t;
 *
 *   int match_function(int v,int h,void* arguments) {
 *     // Extract parameters
 *     match_function_params_t* match_arguments = (match_function_params_t*)arguments;
 *     // Check match
 *     if (v > match_arguments->pattern_length || h > match_arguments->text_length) return 0;
 *     return (match_arguments->pattern[v] == match_arguments->text[h]);
 *   }
 */
typedef int (*alignment_match_funct_t)(int,int,void*);

/*
 * Alignment system configuration
 */
typedef struct {
  // Limits
  int max_alignment_score;       // Maximum score allowed before quit
  // Probing intervals
  int probe_interval_global;     // Score-ticks interval to check any limits
  int probe_interval_compact;    // Score-ticks interval to check BT-buffer compacting
  // Memory
  uint64_t max_partial_compacts; // Maximum partial-compacts before attempting full-compact
  uint64_t max_memory_compact;   // Maximum BT-buffer memory allowed before trigger compact
  uint64_t max_memory_resident;  // Maximum memory allowed to be buffered before reap
  uint64_t max_memory_abort;     // Maximum memory allowed to be used before aborting alignment
  // Verbose
  //  0 - Quiet
  //  1 - Report each sequence aligned                      (brief)
  //  2 - Report each sequence/subsequence aligned          (brief)
  //  3 - Report WFA progress (heavy tasks)                 (verbose)
  //  4 - Full report of each sequence/subsequence aligned  (very verbose)
  int verbose;                   // Verbose (regulates messages during alignment)
  // Debug
  bool check_alignment_correct;  // Verify that the alignment CIGAR output is correct
  // Profile
  profiler_timer_t timer;        // Time alignment
  // OS
  int max_num_threads;           // Maximum number of threads to use to compute/extend WFs
  int min_offsets_per_thread;    // Minimum amount of offsets to spawn a thread
} alignment_system_t;

/*
 * Low-memory modes
 */
typedef enum {
  wavefront_memory_high     = 0, // High-memore mode (fastest, stores all WFs explicitly)
  wavefront_memory_med      = 1, // Succing-memory mode piggyback-based (medium, offloads half-full BT-blocks)
  wavefront_memory_low      = 2, // Succing-memory mode piggyback-based (slow, offloads only full BT-blocks)
  wavefront_memory_ultralow = 3, // Bidirectional WFA
} wavefront_memory_t;

/*
 * Wavefront Aligner Attributes
 */
typedef struct {
  // Distance model
  distance_metric_t distance_metric;       // Alignment metric/distance used
  alignment_scope_t alignment_scope;       // Alignment scope (score only or full-CIGAR)
  alignment_form_t alignment_form;         // Alignment mode (end-to-end/ends-free)
  // Penalties
  linear_penalties_t linear_penalties;     // Gap-linear penalties (placeholder)
  affine_penalties_t affine_penalties;     // Gap-affine penalties (placeholder)
  affine2p_penalties_t affine2p_penalties; // Gap-affine-2p penalties (placeholder)
  // Heuristic strategy
  wavefront_heuristic_t heuristic;         // Wavefront heuristic
  // Memory model
  wavefront_memory_t memory_mode;          // Wavefront memory strategy (modular wavefronts and piggyback)
  // Custom function to compare sequences
  alignment_match_funct_t match_funct;     // Custom matching function (match(v,h,args))
  void* match_funct_arguments;             // Generic arguments passed to matching function (args)
  // External MM (instead of allocating one inside)
  mm_allocator_t* mm_allocator;            // MM-Allocator
  // Display
  wavefront_plot_attr_t plot;              // Plot wavefront
  // System
  alignment_system_t system;               // System related parameters
} wavefront_aligner_attr_t;

/*
 * Default parameters
 */
extern wavefront_aligner_attr_t wavefront_aligner_attr_default;

#endif /* WAVEFRONT_ATTRIBUTES_H_ */
