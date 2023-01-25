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
 * DESCRIPTION: WaveFront alignment module for computing wavefronts (edit/indel)
 */

#include "../utils/string_padded.h"
#include "wavefront_compute.h"
#include "wavefront_backtrace_offload.h"

#ifdef WFA_PARALLEL
#include <omp.h>
#endif

/*
 * Compute Kernels
 */
void wavefront_compute_indel_idm(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wf_prev,
    wavefront_t* const wf_curr,
    const int lo,
    const int hi) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  const wf_offset_t* const prev_offsets = wf_prev->offsets;
  wf_offset_t* const curr_offsets = wf_curr->offsets;
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE
  for (k=lo;k<=hi;++k) {
    // Compute maximum offset
    const wf_offset_t ins = prev_offsets[k-1] + 1;
    const wf_offset_t del = prev_offsets[k+1];
    wf_offset_t max = MAX(del,ins);
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    curr_offsets[k] = max;
  }
}
void wavefront_compute_edit_idm(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wf_prev,
    wavefront_t* const wf_curr,
    const int lo,
    const int hi) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  const wf_offset_t* const prev_offsets = wf_prev->offsets;
  wf_offset_t* const curr_offsets = wf_curr->offsets;
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE
  for (k=lo;k<=hi;++k) {
    // Compute maximum offset
    const wf_offset_t ins = prev_offsets[k-1]; // Lower
    const wf_offset_t del = prev_offsets[k+1]; // Upper
    const wf_offset_t misms = prev_offsets[k]; // Mid
    wf_offset_t max = MAX(del,MAX(ins,misms)+1);
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    curr_offsets[k] = max;
  }
}
/*
 * Compute Kernel (Piggyback)
 */
void wavefront_compute_indel_idm_piggyback(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wf_prev,
    wavefront_t* const wf_curr,
    const int lo,
    const int hi,
    const int score) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // Previous WF
  const wf_offset_t* const prev_offsets = wf_prev->offsets;
  const pcigar_t* const prev_pcigar = wf_prev->bt_pcigar;
  const bt_block_idx_t* const prev_bt_idx = wf_prev->bt_prev;
  // Current WF
  wf_offset_t* const curr_offsets = wf_curr->offsets;
  pcigar_t* const curr_pcigar = wf_curr->bt_pcigar;
  bt_block_idx_t* const curr_bt_idx = wf_curr->bt_prev;
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE // Ifs predicated by the compiler
  for (k=lo;k<=hi;++k) {
    // Compute maximum offset
    const wf_offset_t ins = prev_offsets[k-1] + 1;
    const wf_offset_t del = prev_offsets[k+1];
    wf_offset_t max = MAX(del,ins);
    // Update pcigar & bt-block
    if (max == del) {
      curr_pcigar[k] = PCIGAR_PUSH_BACK_DEL(prev_pcigar[k+1]);
      curr_bt_idx[k] = prev_bt_idx[k+1];
    } else { // max == ins
      curr_pcigar[k] = PCIGAR_PUSH_BACK_INS(prev_pcigar[k-1]);
      curr_bt_idx[k] = prev_bt_idx[k-1];
    }
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    curr_offsets[k] = max;
  }
}
void wavefront_compute_edit_idm_piggyback(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wf_prev,
    wavefront_t* const wf_curr,
    const int lo,
    const int hi,
    const int score) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // Previous WF
  const wf_offset_t* const prev_offsets = wf_prev->offsets;
  const pcigar_t* const prev_pcigar = wf_prev->bt_pcigar;
  const bt_block_idx_t* const prev_bt_idx = wf_prev->bt_prev;
  // Current WF
  wf_offset_t* const curr_offsets = wf_curr->offsets;
  pcigar_t* const curr_pcigar = wf_curr->bt_pcigar;
  bt_block_idx_t* const curr_bt_idx = wf_curr->bt_prev;
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE // Ifs predicated by the compiler
  for (k=lo;k<=hi;++k) {
    // Compute maximum offset
    const wf_offset_t ins = prev_offsets[k-1] + 1; // Lower
    const wf_offset_t del = prev_offsets[k+1];     // Upper
    const wf_offset_t misms = prev_offsets[k] + 1; // Mid
    wf_offset_t max = MAX(del,MAX(ins,misms));
    // Update pcigar & bt-block
    if (max == ins) {
      curr_pcigar[k] = PCIGAR_PUSH_BACK_INS(prev_pcigar[k-1]);
      curr_bt_idx[k] = prev_bt_idx[k-1];
    }
    if (max == del) {
      curr_pcigar[k] = PCIGAR_PUSH_BACK_DEL(prev_pcigar[k+1]);
      curr_bt_idx[k] = prev_bt_idx[k+1];
    }
    if (max == misms) {
      curr_pcigar[k] = PCIGAR_PUSH_BACK_MISMS(prev_pcigar[k]);
      curr_bt_idx[k] = prev_bt_idx[k];
    }
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    curr_offsets[k] = max;
  }
}
/*
 * Exact pruning paths
 */
int wf_compute_edit_best_score(
    const int pattern_length,
    const int text_length,
    const int k,
    const wf_offset_t offset) {
  // Compute best-alignment case
  const int left_v = pattern_length - WAVEFRONT_V(k,offset);
  const int left_h = text_length - WAVEFRONT_H(k,offset);
  return (left_v >= left_h) ? left_v - left_h : left_h - left_v;
}
int wf_compute_edit_worst_score(
    const int pattern_length,
    const int text_length,
    const int k,
    const wf_offset_t offset) {
  // Compute worst-alignment case
  const int left_v = pattern_length - WAVEFRONT_V(k,offset);
  const int left_h = text_length - WAVEFRONT_H(k,offset);
  return MAX(left_v,left_h);
}
void wavefront_compute_edit_exact_prune(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wavefront) {
  // Parameters
  const int plen = wf_aligner->pattern_length;
  const int tlen = wf_aligner->text_length;
  wf_offset_t* const offsets = wavefront->offsets;
  const int lo = wavefront->lo;
  const int hi = wavefront->hi;
  // Speculative compute if needed
  if (WAVEFRONT_LENGTH(lo,hi) < 1000) return;
  const int sample_k = lo + (hi-lo)/2;
  const wf_offset_t sample_offset = offsets[sample_k];
  if (sample_offset < 0) return; // Unlucky null in the middle
  const int smax_sample = wf_compute_edit_worst_score(plen,tlen,sample_k,offsets[sample_k]);
  const int smin_lo = wf_compute_edit_best_score(plen,tlen,lo,offsets[lo]);
  const int smin_hi = wf_compute_edit_best_score(plen,tlen,hi,offsets[hi]);
  if (smin_lo <= smax_sample && smin_hi <= smax_sample) return;
  /*
   * Suggested by Heng Li as an effective exact-prunning technique
   * for sequences of very different length where some diagonals
   * can be proven impossible to yield better alignments.
   */
  // Compute the best worst-case-alignment
  int score_min_worst = INT_MAX;
  int k;
  for (k=lo;k<=hi;++k) {
    const wf_offset_t offset = offsets[k];
    if (offset < 0) continue; // Skip nulls
    // Compute worst-alignment case
    const int score_worst = wf_compute_edit_worst_score(plen,tlen,k,offset);
    if (score_worst < score_min_worst) score_min_worst = score_worst;
  }
  // Compare against the best-case-alignment (Prune from bottom)
  int lo_reduced = lo;
  for (k=lo;k<=hi;++k) {
    // Compute best-alignment case
    const wf_offset_t offset = offsets[k];
    const int score_best = wf_compute_edit_best_score(plen,tlen,k,offset);
    // Compare best and worst
    if (score_best <= score_min_worst) break;
    ++lo_reduced;
  }
  wavefront->lo = lo_reduced;
  // Compare against the best-case-alignment (Prune from top)
  int hi_reduced = hi;
  for (k=hi;k>lo_reduced;--k) {
    // Compute best-alignment case
    const wf_offset_t offset = offsets[k];
    const int score_best = wf_compute_edit_best_score(plen,tlen,k,offset);
    // Compare best and worst
    if (score_best <= score_min_worst) break;
    --hi_reduced;
  }
  wavefront->hi = hi_reduced;
}
/*
 * Compute next wavefront
 */
void wavefront_compute_edit_dispatcher(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    wavefront_t* const wf_prev,
    wavefront_t* const wf_curr,
    const int lo,
    const int hi) {
  if (wf_aligner->wf_components.bt_piggyback) {
    if (wf_aligner->penalties.distance_metric == indel) {
      wavefront_compute_indel_idm_piggyback(wf_aligner,wf_prev,wf_curr,lo,hi,score);
    } else {
      wavefront_compute_edit_idm_piggyback(wf_aligner,wf_prev,wf_curr,lo,hi,score);
    }
  } else {
    if (wf_aligner->penalties.distance_metric == indel) {
      wavefront_compute_indel_idm(wf_aligner,wf_prev,wf_curr,lo,hi);
    } else {
      wavefront_compute_edit_idm(wf_aligner,wf_prev,wf_curr,lo,hi);
    }
  }
}
void wavefront_compute_edit_dispatcher_omp(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const wf_prev,
    wavefront_t* const wf_curr,
    const int lo,
    const int hi,
    const int score) {
  // Parameters
  const int num_threads = wavefront_compute_num_threads(wf_aligner,lo,hi);
  // Multithreading dispatcher
  if (num_threads == 1) {
    // Compute next wavefront
    wavefront_compute_edit_dispatcher(
        wf_aligner,score,wf_prev,wf_curr,lo,hi);
  } else {
#ifdef WFA_PARALLEL
    // Compute next wavefront in parallel
    #pragma omp parallel num_threads(num_threads)
    {
      int t_lo, t_hi;
      const int thread_id = omp_get_thread_num();
      const int thread_num = omp_get_num_threads();
      wavefront_compute_thread_limits(thread_id,thread_num,lo,hi,&t_lo,&t_hi);
      wavefront_compute_edit_dispatcher(
          wf_aligner,score,wf_prev,wf_curr,t_lo,t_hi);
    }
#endif
  }
}
void wavefront_compute_edit(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  // Compute scores
  int score_prev = score - 1;
  int score_curr = score;
  if (wf_components->memory_modular) { // Modular wavefront
    score_prev = score_prev % wf_components->max_score_scope;
    score_curr = score_curr % wf_components->max_score_scope;
    if (wf_components->mwavefronts[score_curr]) { // Free
      wavefront_slab_free(wf_aligner->wavefront_slab,wf_components->mwavefronts[score_curr]);
    }
  }
  // Fetch previous wavefront, compute limits & initialize
  wavefront_t* const wf_prev = wf_components->mwavefronts[score_prev];
  const int lo = wf_prev->lo - 1;
  const int hi = wf_prev->hi + 1;
  //  wf_components->historic_min_lo = min_lo;
  //  wf_components->historic_max_hi = max_hi;
  wf_prev->offsets[lo-1] = WAVEFRONT_OFFSET_NULL;
  wf_prev->offsets[lo] = WAVEFRONT_OFFSET_NULL;
  wf_prev->offsets[hi] = WAVEFRONT_OFFSET_NULL;
  wf_prev->offsets[hi+1] = WAVEFRONT_OFFSET_NULL;
  // Allocate output wavefront
  wavefront_t* const wf_curr = wavefront_slab_allocate(wf_aligner->wavefront_slab,lo-2,hi+2);
  wf_components->mwavefronts[score_curr] = wf_curr;
  wf_components->mwavefronts[score_curr]->lo = lo;
  wf_components->mwavefronts[score_curr]->hi = hi;
  // Compute Wavefront
  wavefront_compute_edit_dispatcher_omp(wf_aligner,wf_prev,wf_curr,lo,hi,score);
  // Offload backtrace (if necessary)
  if (wf_components->bt_piggyback && score % PCIGAR_MAX_LENGTH == 0) {
    wavefront_backtrace_offload_blocks_linear(
        wf_aligner,wf_curr->offsets,wf_curr->bt_pcigar,wf_curr->bt_prev,lo,hi);
  }
  // Trim wavefront ends
  wavefront_compute_trim_ends(wf_aligner,wf_curr);
  if (wf_curr->null) wf_aligner->align_status.num_null_steps = INT_MAX;
  // Exact pruning paths
  if (wf_aligner->alignment_form.span == alignment_end2end &&
      wf_aligner->penalties.distance_metric == edit) {
    wavefront_compute_edit_exact_prune(wf_aligner,wf_curr);
  }
}


