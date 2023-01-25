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
 * DESCRIPTION: Support functions for wavefront heuristic strategies
 */

#include "wavefront_heuristic.h"
#include "wavefront_aligner.h"

/*
 * Setup
 */
void wavefront_heuristic_set_none(
    wavefront_heuristic_t* const wf_heuristic) {
  wf_heuristic->strategy = wf_heuristic_none;
}
void wavefront_heuristic_set_wfadaptive(
    wavefront_heuristic_t* const wf_heuristic,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int steps_between_cutoffs) {
  wf_heuristic->strategy |= wf_heuristic_wfadaptive;
  wf_heuristic->min_wavefront_length = min_wavefront_length;
  wf_heuristic->max_distance_threshold = max_distance_threshold;
  wf_heuristic->steps_between_cutoffs = steps_between_cutoffs;
  // Internals
  wf_heuristic->steps_wait = steps_between_cutoffs;
}
void wavefront_heuristic_set_wfmash(
    wavefront_heuristic_t* const wf_heuristic,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int steps_between_cutoffs) {
  wf_heuristic->strategy |= wf_heuristic_wfmash;
  wf_heuristic->min_wavefront_length = min_wavefront_length;
  wf_heuristic->max_distance_threshold = max_distance_threshold;
  wf_heuristic->steps_between_cutoffs = steps_between_cutoffs;
  // Internals
  wf_heuristic->steps_wait = steps_between_cutoffs;
}
void wavefront_heuristic_set_xdrop(
    wavefront_heuristic_t* const wf_heuristic,
    const int xdrop,
    const int steps_between_cutoffs) {
  wf_heuristic->strategy |= wf_heuristic_xdrop;
  wf_heuristic->xdrop = xdrop;
  wf_heuristic->steps_between_cutoffs = steps_between_cutoffs;
  // Internals
  wf_heuristic->steps_wait = steps_between_cutoffs;
  wf_heuristic->max_sw_score = 0;
  wf_heuristic->max_sw_score_offset = WAVEFRONT_OFFSET_NULL;
  wf_heuristic->max_sw_score_k = DPMATRIX_DIAGONAL_NULL;
}
void wavefront_heuristic_set_zdrop(
    wavefront_heuristic_t* const wf_heuristic,
    const int zdrop,
    const int steps_between_cutoffs) {
  wf_heuristic->strategy |= wf_heuristic_zdrop;
  wf_heuristic->zdrop = zdrop;
  wf_heuristic->steps_between_cutoffs = steps_between_cutoffs;
  // Internals
  wf_heuristic->steps_wait = steps_between_cutoffs;
  wf_heuristic->max_sw_score = 0;
  wf_heuristic->max_sw_score_offset = WAVEFRONT_OFFSET_NULL;
  wf_heuristic->max_sw_score_k = DPMATRIX_DIAGONAL_NULL;
}
void wavefront_heuristic_set_banded_static(
    wavefront_heuristic_t* const wf_heuristic,
    const int band_min_k,
    const int band_max_k) {
  wf_heuristic->strategy |= wf_heuristic_banded_static;
  wf_heuristic->min_k = band_min_k;
  wf_heuristic->max_k = band_max_k;
}
void wavefront_heuristic_set_banded_adaptive(
    wavefront_heuristic_t* const wf_heuristic,
    const int band_min_k,
    const int band_max_k,
    const int steps_between_cutoffs) {
  wf_heuristic->strategy |= wf_heuristic_banded_adaptive;
  wf_heuristic->min_k = band_min_k;
  wf_heuristic->max_k = band_max_k;
  wf_heuristic->steps_between_cutoffs = steps_between_cutoffs;
  // Internals
  wf_heuristic->steps_wait = steps_between_cutoffs;
}
void wavefront_heuristic_clear(
    wavefront_heuristic_t* const wf_heuristic) {
  // Internals
  wf_heuristic->steps_wait = wf_heuristic->steps_between_cutoffs;
  wf_heuristic->max_sw_score = 0;
  wf_heuristic->max_sw_score_offset = WAVEFRONT_OFFSET_NULL;
  wf_heuristic->max_sw_score_k = DPMATRIX_DIAGONAL_NULL;
}
/*
 * Utils
 */
int wf_distance_end2end(
    const wf_offset_t offset,
    const int k,
    const int pattern_length,
    const int text_length) {
  const int left_v = pattern_length - WAVEFRONT_V(k,offset);
  const int left_h = text_length - WAVEFRONT_H(k,offset);
  return (offset >= 0) ? MAX(left_v,left_h) : -WAVEFRONT_OFFSET_NULL;
}
int wf_distance_end2end_weighted(
    const wf_offset_t offset,
    const int k,
    const int pattern_length,
    const int text_length,
    const int mfactor) {
  const int v = WAVEFRONT_V(k,offset);
  const int h = WAVEFRONT_H(k,offset);
  const int left_v = ((float)(pattern_length - v)/pattern_length * mfactor);
  const int left_h = ((float)(text_length - h)/text_length * mfactor);
  return (offset >= 0) ? MAX(left_v,left_h) : -WAVEFRONT_OFFSET_NULL;
}
int wf_distance_endsfree(
    const wf_offset_t offset,
    const int k,
    const int pattern_length,
    const int text_length,
    const int pattern_end_free,
    const int text_end_free) {
  const int left_v = pattern_length - WAVEFRONT_V(k,offset);
  const int left_h = text_length - WAVEFRONT_H(k,offset);
  const int left_v_endsfree = left_v - pattern_end_free;
  const int left_h_endsfree = left_h - text_end_free;
  const int dist_up = MAX(left_h,left_v_endsfree);
  const int dist_down = MAX(left_v,left_h_endsfree);
  return (offset >= 0) ? MIN(dist_up,dist_down) : -WAVEFRONT_OFFSET_NULL;
}
void wf_heuristic_equate(
    wavefront_t* const wavefront_dst,
    wavefront_t* const wavefront_src) {
  if (wavefront_dst != NULL) {
    if (wavefront_src->lo > wavefront_dst->lo) wavefront_dst->lo = wavefront_src->lo;
    if (wavefront_src->hi < wavefront_dst->hi) wavefront_dst->hi = wavefront_src->hi;
    if (wavefront_dst->lo > wavefront_dst->hi) wavefront_dst->null = true;
    // Save min/max WF initialized
    wavefront_dst->wf_elements_init_min = wavefront_dst->lo;
    wavefront_dst->wf_elements_init_max = wavefront_dst->hi;
  }
}
/*
 * Heuristic Cut-off Wavefront-Adaptive
 */
int wf_compute_distance_end2end(
    wavefront_t* const wavefront,
    const int pattern_length,
    const int text_length,
    wf_offset_t* const distances) {
  // Compute min-distance
  const wf_offset_t* const offsets = wavefront->offsets;
  int k, min_distance = MAX(pattern_length,text_length);
  PRAGMA_LOOP_VECTORIZE
  for (k=wavefront->lo;k<=wavefront->hi;++k) {
    const int distance = wf_distance_end2end(
        offsets[k],k,pattern_length,text_length);
    distances[k] = distance;
    min_distance = MIN(min_distance,distance);
  }
  return min_distance;
}
int wf_compute_distance_end2end_weighted(
    wavefront_t* const wavefront,
    const int pattern_length,
    const int text_length,
    wf_offset_t* const distances) {
  // Parameters
  const int mfactor = ((float)(pattern_length + text_length) / 2); // Mean sequence length
  // Compute min-distance
  const wf_offset_t* const offsets = wavefront->offsets;
  int k, min_distance = MAX(pattern_length,text_length);
  PRAGMA_LOOP_VECTORIZE
  for (k=wavefront->lo;k<=wavefront->hi;++k) {
    const int distance = wf_distance_end2end_weighted(
        offsets[k],k,pattern_length,text_length,mfactor);
    distances[k] = distance;
    min_distance = MIN(min_distance,distance);
  }
  return min_distance;
}
int wf_compute_distance_endsfree(
    wavefront_t* const wavefront,
    const int pattern_length,
    const int text_length,
    const int pattern_end_free,
    const int text_end_free,
    wf_offset_t* const distances) {
  // Compute min-distance
  const wf_offset_t* const offsets = wavefront->offsets;
  int k, min_distance = MAX(pattern_length,text_length);
  PRAGMA_LOOP_VECTORIZE
  for (k=wavefront->lo;k<=wavefront->hi;++k) {
    const int distance = wf_distance_endsfree(
        offsets[k],k,pattern_length,text_length,
        pattern_end_free,text_end_free);
    distances[k] = distance;
    min_distance = MIN(min_distance,distance);
  }
  return min_distance;
}
void wf_heuristic_wfadaptive_reduce(
    wavefront_t* const wavefront,
    const wf_offset_t* const distances,
    const int min_distance,
    const int max_distance_threshold,
    const int min_k,
    const int max_k) {
  int k;
  // Reduce from bottom
  const int top_limit = MIN(max_k,wavefront->hi); // Preserve target-diagonals
  int lo_reduced = wavefront->lo;
  for (k=wavefront->lo;k<top_limit;++k) {
    if (distances[k] - min_distance  <= max_distance_threshold) break;
    ++lo_reduced;
  }
  wavefront->lo = lo_reduced;
  // Reduce from top
  const int botton_limit = MAX(min_k,wavefront->lo); // Preserve target-diagonals
  int hi_reduced = wavefront->hi;
  for (k=wavefront->hi;k>botton_limit;--k) {
    if (distances[k] - min_distance <= max_distance_threshold) break;
    --hi_reduced;
  }
  wavefront->hi = hi_reduced;
}
void wavefront_heuristic_wfadaptive(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront,
    const bool wfmash_mode) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  const int min_wavefront_length = wf_aligner->heuristic.min_wavefront_length;
  const int max_distance_threshold = wf_aligner->heuristic.max_distance_threshold;
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  // Check steps
  if (wf_heuristic->steps_wait > 0) return;
  // Check minimum wavefront length
  const int base_hi = wavefront->hi;
  const int base_lo = wavefront->lo;
  if ((base_hi - base_lo + 1) < min_wavefront_length) return;
  // Use victim as temporal buffer
  wavefront_components_resize_null__victim(&wf_aligner->wf_components,base_lo-1,base_hi+1);
  wf_offset_t* const distances = wf_aligner->wf_components.wavefront_victim->offsets;
  // Compute distance & cut-off
  int min_distance;
  if (wfmash_mode) {
    min_distance = wf_compute_distance_end2end_weighted(
        wavefront,pattern_length,text_length,distances);
  } else {
    min_distance = wf_compute_distance_end2end(
        wavefront,pattern_length,text_length,distances);
  }
  // Cut-off wavefront
  const int alignment_k = DPMATRIX_DIAGONAL(text_length,pattern_length);
  wf_heuristic_wfadaptive_reduce(
      wavefront,distances,min_distance,max_distance_threshold,
      alignment_k,alignment_k);
  // Set wait steps (don't repeat this heuristic often)
  wf_heuristic->steps_wait = wf_heuristic->steps_between_cutoffs;
}
/*
 * Heuristic Cut-off Drops
 */
void wf_heuristic_compute_sw_scores(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront,
    const int wf_score,
    wf_offset_t* const sw_scores,
    wf_offset_t* const max_sw_score,
    wf_offset_t* const max_k,
    wf_offset_t* const max_offset) {
  // Parameters
  const int wf_match = wf_aligner->penalties.match;
  const int swg_match = (wf_match==0) ? 1 : -(wf_aligner->penalties.match);
  // Compute min-distance
  const wf_offset_t* const offsets = wavefront->offsets;
  int k, cmax_sw_score = INT_MIN, cmax_k = 0, cmax_offset = 0;
  PRAGMA_LOOP_VECTORIZE
  for (k=wavefront->lo;k<=wavefront->hi;++k) {
    const wf_offset_t offset = offsets[k];
    if (offset < 0) continue;
    const int v = WAVEFRONT_V(k,offset);
    const int h = WAVEFRONT_H(k,offset);
    const int sw_score = (wf_match==0) ?
        (swg_match*(v+h) - wf_score) :
        WF_SCORE_TO_SW_SCORE(swg_match,v,h,wf_score);
    sw_scores[k] = sw_score;
    if (cmax_sw_score < sw_score) {
      cmax_sw_score = sw_score;
      cmax_k = k;
      cmax_offset = offset;
    }
  }
  *max_sw_score = cmax_sw_score;
  *max_k = cmax_k;
  *max_offset = cmax_offset;
}
void wavefront_heuristic_xdrop(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront,
    const int score) {
  // Parameters
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  const int base_hi = wavefront->hi;
  const int base_lo = wavefront->lo;
  // Check steps
  if (wf_heuristic->steps_wait > 0) return;
  // Use victim as temporal buffer
  wavefront_components_resize_null__victim(&wf_aligner->wf_components,base_lo-1,base_hi+1);
  wf_offset_t* const sw_scores = wf_aligner->wf_components.wavefront_victim->offsets;
  // Compute SW scores
  wf_offset_t cmax_sw_score, cmax_k, dummy;
  wf_heuristic_compute_sw_scores(
      wf_aligner,wavefront,score,sw_scores,
      &cmax_sw_score,&cmax_k,&dummy);
  // Apply X-Drop
  const int xdrop = wf_heuristic->xdrop;
  const int max_sw_score = wf_heuristic->max_sw_score;
  const wf_offset_t* const offsets = wavefront->offsets;
  if (wf_heuristic->max_sw_score_k != DPMATRIX_DIAGONAL_NULL) {
    // Reduce from bottom
    int k;
    for (k=wavefront->lo;k<=wavefront->hi;++k) {
      if (offsets[k] < 0) continue;
      //fprintf(stderr,"[XDROP] (max=%d,current=%d) diff=%d leeway=%d\n",
      //    max_sw_score,(int)sw_scores[k],
      //    max_sw_score - (int)sw_scores[k],xdrop);
      if (max_sw_score - (int)sw_scores[k] < xdrop) break;
    }
    wavefront->lo = k;
    // Reduce from top
    for (k=wavefront->hi;k>=wavefront->lo;--k) {
      if (offsets[k] < 0) continue;
      //fprintf(stderr,"[XDROP] (max=%d,current=%d) diff=%d leeway=%d\n",
      //    max_sw_score,(int)sw_scores[k],
      //    max_sw_score - (int)sw_scores[k],xdrop);
      if (max_sw_score - (int)sw_scores[k] < xdrop) break;
    }
    wavefront->hi = k;
    // Update maximum score observed
    if (cmax_sw_score > wf_heuristic->max_sw_score) {
      wf_heuristic->max_sw_score = cmax_sw_score;
      wf_heuristic->max_sw_score_k = cmax_k;
    }
  } else {
    // Update maximum score observed
    wf_heuristic->max_sw_score = cmax_sw_score;
    wf_heuristic->max_sw_score_k = cmax_k;
  }
  // Set wait steps (don't repeat this heuristic often)
  wf_heuristic->steps_wait = wf_heuristic->steps_between_cutoffs;
}
int wf_zdrop_gap_score(
    const int gap_extension_penalty,
    const wf_offset_t offset_1,
    const int k_1,
    const wf_offset_t offset_2,
    const int k_2) {
  int diff_h = WAVEFRONT_H(k_2,offset_2) - WAVEFRONT_H(k_1,offset_1);
  if (diff_h < 0) diff_h = -diff_h;
  int diff_v = WAVEFRONT_V(k_2,offset_2) - WAVEFRONT_V(k_1,offset_1);
  if (diff_v < 0) diff_v = -diff_v;
  const int gap_length = (diff_h >= diff_v) ? diff_h-diff_v : diff_v-diff_h;
  return gap_length * gap_extension_penalty;
}
void wavefront_heuristic_zdrop(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront,
    const int score) {
  // Parameters
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  const int base_hi = wavefront->hi;
  const int base_lo = wavefront->lo;
  // Check steps
  if (wf_heuristic->steps_wait > 0) return;
  // Use victim as temporal buffer
  wavefront_components_resize_null__victim(&wf_aligner->wf_components,base_lo-1,base_hi+1);
  wf_offset_t* const sw_scores = wf_aligner->wf_components.wavefront_victim->offsets;
  // Compute SW scores
  wf_offset_t cmax_sw_score, cmax_k, cmax_offset;
  wf_heuristic_compute_sw_scores(
      wf_aligner,wavefront,score,sw_scores,
      &cmax_sw_score,&cmax_k,&cmax_offset);
  // Apply Z-Drop
  wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  const int gap_e = (penalties->gap_extension1 > 0) ? penalties->gap_extension1 : 1;
  const int zdrop = wf_heuristic->zdrop;
  const int max_sw_score = wf_heuristic->max_sw_score;
  const int max_k = wf_heuristic->max_sw_score_k;
  const int max_offset = wf_heuristic->max_sw_score_offset;
  if (max_k != DPMATRIX_DIAGONAL_NULL) {
    // Update maximum score observed
    if (cmax_sw_score > wf_heuristic->max_sw_score) {
      wf_heuristic->max_sw_score = cmax_sw_score;
      wf_heuristic->max_sw_score_k = cmax_k;
      wf_heuristic->max_sw_score_offset = cmax_offset;
    } else {
      // Test Z-drop
      const int gap_score = wf_zdrop_gap_score(gap_e,max_offset,max_k,cmax_offset,cmax_k);
      //  fprintf(stderr,"[Z-DROP] (max=%d~(%d,%d),current=%d~(%d,%d)) diff=%d leeway=%d\n",
      //      max_sw_score,WAVEFRONT_V(max_k,max_offset),WAVEFRONT_H(max_k,max_offset),
      //      cmax_sw_score,WAVEFRONT_V(cmax_k,cmax_offset),WAVEFRONT_H(cmax_k,cmax_offset),
      //      max_sw_score - cmax_sw_score,
      //      zdrop + gap_score);
      if (max_sw_score - (int)cmax_sw_score > zdrop + gap_score) {
        wavefront->lo = wavefront->hi + 1;
        return; // Z-dropped
      }
    }
  } else {
    // Update maximum score observed
    wf_heuristic->max_sw_score = cmax_sw_score;
    wf_heuristic->max_sw_score_k = cmax_k;
    wf_heuristic->max_sw_score_offset = cmax_offset;
  }
  // Set wait steps (don't repeat this heuristic often)
  wf_heuristic->steps_wait = wf_heuristic->steps_between_cutoffs;
}
/*
 * Heuristic Cut-off Banded
 */
void wavefront_heuristic_banded_static(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront) {
  // Parameters
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  // Check wavefront limits
  if (wavefront->lo < wf_heuristic->min_k) wavefront->lo = wf_heuristic->min_k;
  if (wavefront->hi > wf_heuristic->max_k) wavefront->hi = wf_heuristic->max_k;
}
void wavefront_heuristic_banded_adaptive(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  // Check steps
  if (wf_heuristic->steps_wait > 0) return;
  // Check wavefront length
  const int lo = wavefront->lo;
  const int hi = wavefront->hi;
  const int wf_length = hi - lo + 1;
  if (wf_length < 4) return; // We cannot do anything here
  // Adjust the band
  const wf_offset_t* const offsets = wavefront->offsets;
  const int max_wf_length = wf_heuristic->max_k - wf_heuristic->min_k + 1;
  if (wf_length > max_wf_length) {
    // Sample wavefront
    const int leeway = (wf_length - max_wf_length) / 2;
    const int quarter = wf_length / 4;
    const int dist_p0 = wf_distance_end2end(
        offsets[lo],lo,pattern_length,text_length);
    const int dist_p1 = wf_distance_end2end(
        offsets[lo+quarter],lo+quarter,pattern_length,text_length);
    const int dist_p2 = wf_distance_end2end(
        offsets[lo+2*quarter],lo+2*quarter,pattern_length,text_length);
    const int dist_p3 = wf_distance_end2end(
        offsets[hi],hi,pattern_length,text_length);
    // Heuristically decide where to place the band
    int new_lo = lo;
    if (dist_p0 > dist_p3) new_lo += leeway;
    if (dist_p1 > dist_p2) new_lo += leeway;
    // Set wavefront limits
    wavefront->lo = new_lo;
    if (wavefront->lo < lo) wavefront->lo = lo;
    wavefront->hi = new_lo + max_wf_length - 1;
    if (wavefront->hi > hi) wavefront->hi = hi;
  }
  // Set wait steps (don't repeat this heuristic often)
  wf_heuristic->steps_wait = wf_heuristic->steps_between_cutoffs;
}
/*
 * Heuristic Cut-offs dispatcher
 */
void wavefront_heuristic_cufoff(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int score_mod) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  wavefront_heuristic_t* const wf_heuristic = &wf_aligner->heuristic;
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_components->mwavefronts[score_mod];
  if (mwavefront == NULL || mwavefront->lo > mwavefront->hi) return;
  // Decrease wait steps
  --(wf_heuristic->steps_wait);
  // Select heuristic (WF-Adaptive)
  if (wf_heuristic->strategy & wf_heuristic_wfadaptive) {
    wavefront_heuristic_wfadaptive(wf_aligner,mwavefront,false);
  } else if (wf_heuristic->strategy & wf_heuristic_wfmash) {
    wavefront_heuristic_wfadaptive(wf_aligner,mwavefront,true);
  }
  // Select heuristic (Drops)
  if (wf_heuristic->strategy & wf_heuristic_xdrop) {
    wavefront_heuristic_xdrop(wf_aligner,mwavefront,score);
  } else if (wf_heuristic->strategy & wf_heuristic_zdrop) {
    wavefront_heuristic_zdrop(wf_aligner,mwavefront,score);
  }
  // Select heuristic (Banded)
  if (wf_heuristic->strategy & wf_heuristic_banded_static) {
    wavefront_heuristic_banded_static(wf_aligner,mwavefront);
  } else if (wf_heuristic->strategy & wf_heuristic_banded_adaptive) {
    wavefront_heuristic_banded_adaptive(wf_aligner,mwavefront);
  }
  // Check wavefront length
  if (mwavefront->lo > mwavefront->hi) mwavefront->null = true;
  // DEBUG
  // const int wf_length_base = hi_base-lo_base+1;
  // const int wf_length_reduced = mwavefront->hi-mwavefront->lo+1;
  // fprintf(stderr,"[WFA::Heuristic] Heuristic from %d to %d offsets (%2.2f%%)\n",
  //    wf_length_base,wf_length_reduced,100.0f*(float)wf_length_reduced/(float)wf_length_base);
  // Save min/max WF initialized
  mwavefront->wf_elements_init_min = mwavefront->lo;
  mwavefront->wf_elements_init_max = mwavefront->hi;
  // Equate other wavefronts
  if (distance_metric <= gap_linear) return;
  // Cut-off the other wavefronts (same dimensions as M)
  wavefront_t* const i1wavefront = wf_components->i1wavefronts[score_mod];
  wavefront_t* const d1wavefront = wf_components->d1wavefronts[score_mod];
  wf_heuristic_equate(i1wavefront,mwavefront);
  wf_heuristic_equate(d1wavefront,mwavefront);
  if (distance_metric == gap_affine) return;
  wavefront_t* const i2wavefront = wf_components->i2wavefronts[score_mod];
  wavefront_t* const d2wavefront = wf_components->d2wavefronts[score_mod];
  wf_heuristic_equate(i2wavefront,mwavefront);
  wf_heuristic_equate(d2wavefront,mwavefront);
}
/*
 * Display
 */
void wavefront_heuristic_print(
    FILE* const stream,
    wavefront_heuristic_t* const wf_heuristic) {
  // Select heuristic strategy
  if (wf_heuristic->strategy == wf_heuristic_none) {
    fprintf(stream,"(none)");
  } else {
    // WF-Adaptive
    if (wf_heuristic->strategy & wf_heuristic_wfadaptive) {
      fprintf(stream,"(wfadapt,%d,%d,%d)",
          wf_heuristic->min_wavefront_length,
          wf_heuristic->max_distance_threshold,
          wf_heuristic->steps_between_cutoffs);
    } else if (wf_heuristic->strategy & wf_heuristic_wfmash) {
      fprintf(stream,"(wfmash,%d,%d,%d)",
          wf_heuristic->min_wavefront_length,
          wf_heuristic->max_distance_threshold,
          wf_heuristic->steps_between_cutoffs);
    }
    // Drops
    if (wf_heuristic->strategy & wf_heuristic_xdrop) {
      fprintf(stream,"(xdrop,%d,%d)",
          wf_heuristic->xdrop,
          wf_heuristic->steps_between_cutoffs);
    }
    if (wf_heuristic->strategy & wf_heuristic_zdrop) {
      fprintf(stream,"(zdrop,%d,%d)",
          wf_heuristic->zdrop,
          wf_heuristic->steps_between_cutoffs);
    }
    // Banded
    if (wf_heuristic->strategy & wf_heuristic_banded_static) {
      fprintf(stream,"(banded-static,%d,%d)",
          wf_heuristic->min_k,
          wf_heuristic->max_k);
    }
    if (wf_heuristic->strategy & wf_heuristic_banded_adaptive) {
      fprintf(stream,"(banded-adapt,%d,%d,%d)",
          wf_heuristic->min_k,
          wf_heuristic->max_k,
          wf_heuristic->steps_between_cutoffs);
    }
  }
}
