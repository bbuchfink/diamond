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

#include "wavefront_attributes.h"

/*
 * Default parameters
 */
wavefront_aligner_attr_t wavefront_aligner_attr_default = {
    // Distance model & Penalties
    .distance_metric = gap_affine,
    .alignment_scope = compute_alignment,
    .alignment_form = {
        .span = alignment_end2end,
        .pattern_begin_free = 0,
        .pattern_end_free = 0,
        .text_begin_free = 0,
        .text_end_free = 0,
    },
    // Custom matching functions
    .match_funct = NULL,           // Use default match-compare function
    .match_funct_arguments = NULL, // No arguments
    // Penalties
    .linear_penalties = {
        .match = 0,
        .mismatch = 4,
        .indel = 2,
    },
    .affine_penalties = {
        .match = 0,
        .mismatch = 4,
        .gap_opening = 6,
        .gap_extension = 2,
    },
    .affine2p_penalties = {
        .match = 0,
        .mismatch = 4,
        .gap_opening1 = 6,
        .gap_extension1 = 2,
        .gap_opening2 = 24,
        .gap_extension2 = 1,
    },
    // Heuristic
    .heuristic = {
        .strategy = wf_heuristic_wfadaptive,
        .min_wavefront_length = 10,
        .max_distance_threshold = 50,
        .steps_between_cutoffs = 1,
    },
    // Memory model
    .memory_mode = wavefront_memory_high,
    // MM
    .mm_allocator = NULL, // Use private MM
    // Display
    .plot = {
        .enabled = false,
        .resolution_points = 2000,
        .align_level = 0,
    },
    // System
    .system = {
        .max_alignment_score = INT_MAX, // Unlimited
        .probe_interval_global = 3000,
        .probe_interval_compact = 6000,
        .max_memory_compact = -1,  // Automatically set based on memory-mode
        .max_memory_resident = -1, // Automatically set based on memory-mode
        .max_memory_abort = UINT64_MAX, // Unlimited
        .verbose = 0, // Quiet
        .check_alignment_correct = false,
        .max_num_threads = 1,           // Single thread by default
        .min_offsets_per_thread = 500   // Minimum WF-length to spawn a thread
    },
};
