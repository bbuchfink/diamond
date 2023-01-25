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

#include "wavefront_plot.h"
#include "wavefront_aligner.h"

/*
 * Heatmaps
 */
void wavefront_plot_heatmaps_allocate(
    wavefront_plot_t* const wf_plot,
    const int pattern_length,
    const int text_length) {
  wavefront_plot_attr_t* const attributes = &wf_plot->attributes;
  // Compute dimensions
  const int resolution_points = attributes->resolution_points;
  const int min_v = (wf_plot->min_v == -1) ? 0 : wf_plot->min_v;
  const int max_v = (wf_plot->max_v == -1) ? pattern_length-1 : wf_plot->max_v;
  const int min_h = (wf_plot->min_h == -1) ? 0 : wf_plot->min_h;
  const int max_h = (wf_plot->max_h == -1) ? text_length-1 : wf_plot->max_h;
  // Behavior
  wf_plot->behavior_heatmap = heatmap_new(heatmap_value,
      min_v,max_v,min_h,max_h,resolution_points);
  // Wavefront Components
  wf_plot->m_heatmap = heatmap_new(heatmap_min,
      min_v,max_v,min_h,max_h,resolution_points);
  wf_plot->i1_heatmap = NULL;
  wf_plot->d1_heatmap = NULL;
  wf_plot->i2_heatmap = NULL;
  wf_plot->d2_heatmap = NULL;
  if (wf_plot->distance_metric < gap_affine) return;
  // Gap-affine
  wf_plot->i1_heatmap = heatmap_new(heatmap_min,
      min_v,max_v,min_h,max_h,resolution_points);
  wf_plot->d1_heatmap = heatmap_new(heatmap_min,
      min_v,max_v,min_h,max_h,resolution_points);
  if (wf_plot->distance_metric == gap_affine) return;
  // Gap-affine-2p
  wf_plot->i2_heatmap = heatmap_new(heatmap_min,
      min_v,max_v,min_h,max_h,resolution_points);
  wf_plot->d2_heatmap = heatmap_new(heatmap_min,
      min_v,max_v,min_h,max_h,resolution_points);
}
void wavefront_plot_heatmaps_free(
    wavefront_plot_t* const wf_plot) {
  heatmap_delete(wf_plot->behavior_heatmap);
  heatmap_delete(wf_plot->m_heatmap);
  if (wf_plot->i1_heatmap) heatmap_delete(wf_plot->i1_heatmap);
  if (wf_plot->d1_heatmap) heatmap_delete(wf_plot->d1_heatmap);
  if (wf_plot->i2_heatmap) heatmap_delete(wf_plot->i2_heatmap);
  if (wf_plot->d2_heatmap) heatmap_delete(wf_plot->d2_heatmap);
}
/*
 * Setup
 */
wavefront_plot_t* wavefront_plot_new(
    const distance_metric_t distance_metric,
    const int pattern_length,
    const int text_length,
    wavefront_plot_attr_t* const attributes) {
  // Handler
  wavefront_plot_t* const wf_plot = (wavefront_plot_t*)malloc(sizeof(wavefront_plot_t));
  // Parameters
  wf_plot->attributes = *attributes;
  wf_plot->distance_metric = distance_metric;
  wf_plot->min_v = -1;
  wf_plot->max_v = -1;
  wf_plot->min_h = -1;
  wf_plot->max_h = -1;
  // Allocate and configure
  wavefront_plot_heatmaps_allocate(wf_plot,pattern_length,text_length);
  // Clear offsets
  wf_plot->offset_h = 0;
  wf_plot->offset_v = 0;
  // Return
  return wf_plot;
}
void wavefront_plot_resize(
    wavefront_plot_t* const wf_plot,
    const int pattern_length,
    const int text_length) {
  // Free heatmaps
  wavefront_plot_heatmaps_free(wf_plot);
  // Allocate new heatmaps
  wavefront_plot_heatmaps_allocate(wf_plot,pattern_length,text_length);
  // Clear offsets
  wf_plot->offset_h = 0;
  wf_plot->offset_v = 0;
}
void wavefront_plot_delete(
    wavefront_plot_t* const wf_plot) {
  // Heatmaps
  wavefront_plot_heatmaps_free(wf_plot);
  // Handler
  free(wf_plot);
}
/*
 * Accessors
 */
void wavefront_plot_component(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront,
    const int score,
    heatmap_t* const wf_heatmap,
    const bool extend) {
  // Check wavefront
  if (wavefront == NULL) return;
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  const char* const pattern = wf_aligner->pattern;
  const char* const text = wf_aligner->text;
  wavefront_plot_t* const plot = wf_aligner->plot;
  const bool reverse = (wf_aligner->align_mode == wf_align_biwfa_breakpoint_reverse);
  // Traverse all offsets
  int k;
  for (k=wavefront->lo;k<=wavefront->hi;++k) {
    const wf_offset_t offset = wavefront->offsets[k];
    if (offset < 0) continue;
    // Compute local coordinates
    int v_local = WAVEFRONT_V(k,offset);
    int h_local = WAVEFRONT_H(k,offset);
    if (v_local < 0 || v_local >= pattern_length) continue;
    if (h_local < 0 || h_local >= text_length) continue;
    // Compute global coordinates
    int v_global, h_global;
    if (reverse) {
      v_global = plot->offset_v + (pattern_length - 1 - v_local);
      h_global = plot->offset_h + (text_length - 1 - h_local);
    } else {
      v_global = plot->offset_v + v_local;
      h_global = plot->offset_h + h_local;
    }
    // Plot
    if (reverse) {
      if (h_local>0 && v_local>0) heatmap_set(wf_heatmap,v_global+1,h_global+1,score);
    } else {
      if (h_local>0 && v_local>0) heatmap_set(wf_heatmap,v_global-1,h_global-1,score);
    }
    // Simulate extension
    if (extend) {
      while (v_local < pattern_length &&
             h_local < text_length &&
             pattern[v_local] == text[h_local]) {
        if (reverse) {
          v_global--; h_global--;
        } else {
          v_global++; h_global++;
        }
        v_local++; h_local++;
        if (reverse) {
          heatmap_set(wf_heatmap,v_global+1,h_global+1,score);
        } else {
          heatmap_set(wf_heatmap,v_global-1,h_global-1,score);
        }
      }
    }
  }
}
void wavefront_plot(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int align_level) {
  //if (wf_aligner->align_mode == wf_align_biwfa_breakpoint_forward) return;
  // Check plotting enabled wrt align-level
  if (wf_aligner->align_mode == wf_align_biwfa_breakpoint_forward ||
      wf_aligner->align_mode == wf_align_biwfa_breakpoint_reverse) {
    if (align_level != wf_aligner->plot->attributes.align_level) return;
  }
  if (wf_aligner->align_mode == wf_align_biwfa_subsidiary &&
      wf_aligner->plot->attributes.align_level != -1) return;
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const int score_mod = (wf_components->memory_modular) ? score%wf_components->max_score_scope : score;
  // Plot wavefront components
  wavefront_plot_component(wf_aligner,
      wf_components->mwavefronts[score_mod],
      score,wf_aligner->plot->m_heatmap,true);
  if (distance_metric < gap_affine) return;
  // Gap-affine
  wavefront_plot_component(wf_aligner,
      wf_components->i1wavefronts[score_mod],
      score,wf_aligner->plot->i1_heatmap,false);
  wavefront_plot_component(wf_aligner,
      wf_components->d1wavefronts[score_mod],
      score,wf_aligner->plot->d1_heatmap,false);
  if (distance_metric == gap_affine) return;
  // Gap-affine-2p
  wavefront_plot_component(wf_aligner,
      wf_components->i2wavefronts[score_mod],
      score,wf_aligner->plot->i2_heatmap,false);
  wavefront_plot_component(wf_aligner,
      wf_components->d2wavefronts[score_mod],
      score,wf_aligner->plot->d2_heatmap,false);
}
/*
 * Display
 */
void wavefront_plot_print_cigar(
    FILE* const stream,
    cigar_t* const cigar,
    const char target_operation) {
  int i, h=0, v=0, count=0;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    // Check operation
    const char operation = cigar->operations[i];
    switch (operation) {
      case 'M': case 'X': ++h; ++v; break;
      case 'I': ++h; break;
      case 'D': ++v; break;
      default: break;
    }
    // Print point
    if (operation == target_operation && h>0 && v>0) {
      if (count++ > 0) fprintf(stream,";");
      fprintf(stream,"%d,%d",h-1,v-1);
    }
  }
}
void wavefront_plot_print(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  wavefront_plot_t* const wf_plot = wf_aligner->plot;
  // Metadata
  if (wf_aligner->match_funct != NULL) {
    fprintf(stream,"# PatternLength 0\n");
    fprintf(stream,"# TextLength 0\n");
    fprintf(stream,"# Pattern -\n");
    fprintf(stream,"# Text -\n");
  } else {
    fprintf(stream,"# PatternLength %d\n",wf_aligner->pattern_length);
    fprintf(stream,"# Pattern %.*s\n",wf_aligner->pattern_length,wf_aligner->pattern);
    fprintf(stream,"# TextLength %d\n",wf_aligner->text_length);
    fprintf(stream,"# Text %.*s\n",wf_aligner->text_length,wf_aligner->text);
  }
  fprintf(stream,"# Penalties ");
  wavefront_penalties_print(stream,&wf_aligner->penalties);
  fprintf(stream,"\n");
  // Alignment mode
  fprintf(stream,"# WFAMode ");
  wavefront_aligner_print_mode(stream,wf_aligner);
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  if (wf_heuristic->strategy != wf_heuristic_none) {
    wavefront_heuristic_print(stream,wf_heuristic);
  }
  fprintf(stream,"\n");
  // Wavefront components
  fprintf(stream,"# Heatmap M\n"); heatmap_print(stream,wf_plot->m_heatmap);
  if (distance_metric == gap_affine) {
    fprintf(stream,"# Heatmap I1\n"); heatmap_print(stream,wf_plot->i1_heatmap);
    fprintf(stream,"# Heatmap D1\n"); heatmap_print(stream,wf_plot->d1_heatmap);
  }
  if (distance_metric == gap_affine_2p) {
    fprintf(stream,"# Heatmap I2\n"); heatmap_print(stream,wf_plot->i2_heatmap);
    fprintf(stream,"# Heatmap D2\n"); heatmap_print(stream,wf_plot->d2_heatmap);
  }
  // Extend
  fprintf(stream,"# Heatmap Extend\n"); heatmap_print(stream,wf_plot->behavior_heatmap);
  // CIGAR
  if (wf_aligner->alignment_scope == compute_alignment) {
    fprintf(stream,"# List CIGAR-M ");
    wavefront_plot_print_cigar(stream,wf_aligner->cigar,'M');
    fprintf(stream,"\n");
    fprintf(stream,"# List CIGAR-X ");
    wavefront_plot_print_cigar(stream,wf_aligner->cigar,'X');
    fprintf(stream,"\n");
    fprintf(stream,"# List CIGAR-I ");
    wavefront_plot_print_cigar(stream,wf_aligner->cigar,'I');
    fprintf(stream,"\n");
    fprintf(stream,"# List CIGAR-D ");
    wavefront_plot_print_cigar(stream,wf_aligner->cigar,'D');
    fprintf(stream,"\n");
  }
}
