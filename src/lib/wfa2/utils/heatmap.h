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
 */

#ifndef HEATMAP_H_
#define HEATMAP_H_

#include "../utils/commons.h"

/*
 * Heatmap
 */
typedef enum {
  heatmap_min,   // Min value stays
  heatmap_max,   // Max value stays
  heatmap_value, // Last value set stays
} heatmap_type;
typedef struct {
  // Configuration
  heatmap_type type;
  // Dimensions
  int num_rows;
  int num_columns;
  // Range
  int min_v;
  int max_v;
  int min_h;
  int max_h;
  float binning_factor;
  // Data
  int** values;
} heatmap_t;

/*
 * Setup
 */
heatmap_t* heatmap_new(
    const heatmap_type type,
    const int min_v,
    const int max_v,
    const int min_h,
    const int max_h,
    const int resolution_points);
void heatmap_clear(
    heatmap_t* const heatmap);
void heatmap_delete(
    heatmap_t* const heatmap);

/*
 * Accessors
 */
void heatmap_set(
    heatmap_t* const heatmap,
    const int v,
    const int h,
    const int value);

/*
 * Display
 */
void heatmap_print(
    FILE* const stream,
    heatmap_t* const heatmap);

#endif /* HEATMAP_H_ */
