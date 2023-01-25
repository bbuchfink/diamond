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

#include "heatmap.h"

/*
 * Constants
 */
#define HEATMAP_INT_MIN INT_MIN
#define HEATMAP_INT_MAX INT_MAX
#define HEATMAP_SEPARATOR ','

/*
 * Setup
 */
heatmap_t* heatmap_new(
    const heatmap_type type,
    const int min_v,
    const int max_v,
    const int min_h,
    const int max_h,
    const int resolution_points) {
  // Alloc
  heatmap_t* const heatmap = (heatmap_t*)malloc(sizeof(heatmap_t));
  // Configure
  heatmap->type = type;
  heatmap->min_v = min_v;
  heatmap->max_v = max_v;
  heatmap->min_h = min_h;
  heatmap->max_h = max_h;
  // Binning (mantain aspect ratio)
  const int v_range = (max_v-min_v+1);
  const int h_range = (max_h-min_h+1);
  const int max_range = MAX(v_range,h_range);
  if (max_range <= resolution_points) {
    heatmap->binning_factor = 1.0f;
    heatmap->num_rows = v_range;
    heatmap->num_columns = h_range;
  } else {
    heatmap->binning_factor = (float)max_range/(float)resolution_points;
    heatmap->num_rows = (float)v_range/heatmap->binning_factor;
    heatmap->num_columns = (float)h_range/heatmap->binning_factor;
  }
  // Allocate matrix
  heatmap->values = (int**)malloc(heatmap->num_rows*sizeof(int*));
  int i;
  for (i=0;i<heatmap->num_rows;++i) {
    heatmap->values[i] = (int*)malloc(heatmap->num_columns*sizeof(int)); // Allocate row
  }
  // Clear
  heatmap_clear(heatmap);
  // Return
  return heatmap;
}
void heatmap_clear(
    heatmap_t* const heatmap) {
  // Parameters
  const heatmap_type type = heatmap->type;
  const int num_rows = heatmap->num_rows;
  const int num_columns = heatmap->num_columns;
  // Clear values
  int i, j;
  for (i=0;i<num_rows;++i) {
    switch (type) {
      case heatmap_min: {
        for (j=0;j<num_columns;++j) heatmap->values[i][j] = HEATMAP_INT_MAX;
        break;
      }
      case heatmap_max:
      case heatmap_value:
      default: {
        for (j=0;j<num_columns;++j) heatmap->values[i][j] = HEATMAP_INT_MIN;
        break;
      }
    }
  }
}
void heatmap_delete(
    heatmap_t* const heatmap) {
  int i;
  for (i=0;i<heatmap->num_rows;++i) {
    free(heatmap->values[i]);
  }
  free(heatmap->values);
  free(heatmap);
}
/*
 * Accessors
 */
void heatmap_set(
    heatmap_t* const heatmap,
    const int v,
    const int h,
    const int value) {
  // Paramters
  const int num_rows = heatmap->num_rows;
  const int num_columns = heatmap->num_columns;
  // Check position
  if (v<heatmap->min_v || v>heatmap->max_v) return;
  if (h<heatmap->min_h || h>heatmap->max_h) return;
  // Adjust position into binning
  int v_adjusted = (float)(v - heatmap->min_v) / heatmap->binning_factor;
  int h_adjusted = (float)(h - heatmap->min_h) / heatmap->binning_factor;
  if (v_adjusted >= num_rows) v_adjusted = num_rows-1;
  if (h_adjusted >= num_columns) h_adjusted = num_columns-1;
  // Set value
  switch (heatmap->type) {
    case heatmap_min:
      heatmap->values[v_adjusted][h_adjusted] = MIN(heatmap->values[v_adjusted][h_adjusted],value);
      break;
    case heatmap_max:
      heatmap->values[v_adjusted][h_adjusted] = MAX(heatmap->values[v_adjusted][h_adjusted],value);
      break;
    case heatmap_value:
    default:
      heatmap->values[v_adjusted][h_adjusted] = value;
      break;
  }
}
/*
 * Display
 */
void heatmap_print(
    FILE* const stream,
    heatmap_t* const heatmap) {
  // Parameters
  const int num_rows = heatmap->num_rows;
  const int num_columns = heatmap->num_columns;
  // Print values
  int v, h;
  for (v=0;v<num_rows;++v) {
    for (h=0;h<num_columns;++h) {
      if (h>0) fprintf(stream,"%c",HEATMAP_SEPARATOR);
      const int value = heatmap->values[v][h];
      if (value == HEATMAP_INT_MIN || value == HEATMAP_INT_MAX) {
        fprintf(stream,"-1");
      } else {
        fprintf(stream,"%d",value);
      }
    }
    fprintf(stream,"\n");
  }
}



