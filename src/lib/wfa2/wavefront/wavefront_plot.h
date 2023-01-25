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
 * DESCRIPTION: Wavefront alignment module for plot
 */

#ifndef WAVEFRONT_PLOT_H_
#define WAVEFRONT_PLOT_H_

#include "../utils/commons.h"
#include "../utils/heatmap.h"
#include "../alignment/score_matrix.h"
#include "../wavefront/wavefront_penalties.h"

// Wavefront ahead definition
typedef struct _wavefront_aligner_t wavefront_aligner_t;

/*
 * Wavefront Display
 */
typedef struct {
  bool enabled;               // Is plotting enabled
  int resolution_points;      // Total resolution points
  int align_level;            // Level of recursion to plot (-1 == final)
} wavefront_plot_attr_t;
typedef struct {
  // Configuration
  wavefront_plot_attr_t attributes;
  distance_metric_t distance_metric;
  int min_v;
  int max_v;
  int min_h;
  int max_h;
  // Wavefront Heatmaps
  heatmap_t* m_heatmap;
  heatmap_t* i1_heatmap;
  heatmap_t* d1_heatmap;
  heatmap_t* i2_heatmap;
  heatmap_t* d2_heatmap;
  heatmap_t* behavior_heatmap;
  // Offsets
  int offset_h;
  int offset_v;
} wavefront_plot_t;

/*
 * Setup
 */
wavefront_plot_t* wavefront_plot_new(
    const distance_metric_t distance_metric,
    const int pattern_length,
    const int text_length,
    wavefront_plot_attr_t* const attributes);
void wavefront_plot_resize(
    wavefront_plot_t* const wf_plot,
    const int pattern_length,
    const int text_length);
void wavefront_plot_delete(
    wavefront_plot_t* const wf_plot);

/*
 * Accessors
 */
void wavefront_plot(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int align_level);

/*
 * Display
 */
void wavefront_plot_print(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner);

#endif /* WAVEFRONT_PLOT_H_ */
