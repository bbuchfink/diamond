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

#pragma once

#include "utils/heatmap.h"
#include "system/profiler_counter.h"
#include "system/profiler_timer.h"
#include "system/mm_stack.h"
#include "alignment/cigar.h"
#include "wfa.h"

/*
 * Initialize wf-alignment conditions
 */
void wavefront_aligner_init_wf(
    wavefront_aligner_t* const wf_aligner);

/*
 * Utils
 */
uint64_t wavefront_aligner_get_size(
    wavefront_aligner_t* const wf_aligner);

void wavefront_aligner_maxtrim_cigar(
    wavefront_aligner_t* const wf_aligner);

/*
 * Display
 */
void wavefront_aligner_print_mode(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_print_scope(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_print_conf(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner);


