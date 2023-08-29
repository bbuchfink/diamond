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
 * DESCRIPTION: WFA module for the "extension" of exact matches
 */

#include "utils/commons.h"
#include "system/mm_allocator.h"
#include "wavefront_extend.h"
#include "wavefront_extend_kernels.h"
#include "wavefront_compute.h"
#include "wavefront_termination.h"

#ifdef WFA_PARALLEL
#include <omp.h>
#endif

/*
 * Wavefront Extension (End-to-end)
 */
void wavefront_extend_end2end_dispatcher_seq(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int lo,
    const int hi) {
  // Parameters
  wavefront_sequences_t* const seqs = &wf_aligner->sequences;
  // Check the sequence mode
  if (seqs->mode == wf_sequences_ascii) {
    wavefront_extend_matches_packed_end2end(wf_aligner,mwavefront,lo,hi);
  } else {
    wf_offset_t dummy;
    wavefront_extend_matches_custom(wf_aligner,mwavefront,score,lo,hi,false,&dummy);
  }
}
void wavefront_extend_end2end_dispatcher_threads(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score) {
  // Parameters
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
  const int num_threads = wavefront_compute_num_threads(wf_aligner,lo,hi);
  if (num_threads == 1) {
    // Extend wavefront single-thread
    wavefront_extend_end2end_dispatcher_seq(wf_aligner,mwavefront,score,lo,hi);
  } else {
#ifdef WFA_PARALLEL
    // Extend wavefront in parallel
    #pragma omp parallel num_threads(num_threads)
    {
      int t_lo, t_hi;
      wavefront_compute_thread_limits(omp_get_thread_num(),omp_get_num_threads(),lo,hi,&t_lo,&t_hi);
      wavefront_extend_end2end_dispatcher_seq(wf_aligner,mwavefront,score,t_lo,t_hi);
    }
#endif
  }
}
int wavefront_extend_end2end(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Compute score
  const bool memory_modular = wf_aligner->wf_components.memory_modular;
  const int max_score_scope = wf_aligner->wf_components.max_score_scope;
  const int score_mod = (memory_modular) ? score % max_score_scope : score;
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score_mod];
  if (mwavefront == NULL) {
    // Check alignment feasibility (for heuristic variants that can lead to no solution)
    if (wf_aligner->align_status.num_null_steps > wf_aligner->wf_components.max_score_scope) {
      wf_aligner->align_status.status = WF_STATUS_UNFEASIBLE;
      wf_aligner->align_status.score = score;
      return 1; // Done
    }
    return 0; // Not done
  }
  // Extend (dispatcher)
  wavefront_extend_end2end_dispatcher_threads(wf_aligner,mwavefront,score);
  const bool end_reached = wavefront_termination_end2end(wf_aligner,mwavefront,score,score_mod);
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_END_REACHED;
    wf_aligner->align_status.score = score;
    return 1; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    if (wavefront_heuristic_cufoff(wf_aligner,score,score_mod)) {
      wf_aligner->align_status.status = WF_STATUS_END_REACHED;
      return 1; // Done
    }
  }
  return 0; // Not done
}
/*
 * Wavefront Extension (End-to-end + MAX-antidiagonal)
 */
wf_offset_t wavefront_extend_end2end_max_dispatcher_seq(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int lo,
    const int hi) {
  // Parameters
  wavefront_sequences_t* const seqs = &wf_aligner->sequences;
  // Check the sequence mode
  if (seqs->mode == wf_sequences_ascii) {
    return wavefront_extend_matches_packed_end2end_max(wf_aligner,mwavefront,lo,hi);
  } else {
    wf_offset_t max_antidiag;
    wavefront_extend_matches_custom(wf_aligner,mwavefront,score,lo,hi,false,&max_antidiag);
    return max_antidiag;
  }
}
wf_offset_t wavefront_extend_end2end_max_dispatcher_threads(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score) {
  // Parameters
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
  wf_offset_t max_antidiag = 0;
  // Select number of threads
  const int num_threads = wavefront_compute_num_threads(wf_aligner,lo,hi);
  if (num_threads == 1) {
    // Extend wavefront single-thread
    max_antidiag = wavefront_extend_end2end_max_dispatcher_seq(wf_aligner,mwavefront,score,lo,hi);
  } else {
    // Extend wavefront in parallel
#ifdef WFA_PARALLEL
    #pragma omp parallel num_threads(num_threads)
    {
      int t_lo, t_hi;
      wavefront_compute_thread_limits(omp_get_thread_num(),omp_get_num_threads(),lo,hi,&t_lo,&t_hi);
      wf_offset_t t_max_antidiag = wavefront_extend_end2end_max_dispatcher_seq(wf_aligner,mwavefront,score,t_lo,t_hi);
      #pragma omp critical
      {
        if (t_max_antidiag > max_antidiag) max_antidiag = t_max_antidiag;
      }
    }
#endif
  }
  // Return maximum antidiagonal
  return max_antidiag;
}
int wavefront_extend_end2end_max(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    int* const max_antidiagonal) {
  // Compute score
  const bool memory_modular = wf_aligner->wf_components.memory_modular;
  const int max_score_scope = wf_aligner->wf_components.max_score_scope;
  const int score_mod = (memory_modular) ? score % max_score_scope : score;
  *max_antidiagonal = 0; // Init
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score_mod];
  if (mwavefront == NULL) {
    // Check alignment feasibility (heuristic variants that can lead to no solution)
    if (wf_aligner->align_status.num_null_steps > wf_aligner->wf_components.max_score_scope) {
      wf_aligner->align_status.status = WF_STATUS_UNFEASIBLE;
      wf_aligner->align_status.score = score;
      return 1; // Done
    }
    return 0; // Not done
  }
  // Extend (dispatcher)
  const wf_offset_t max_ak = wavefront_extend_end2end_max_dispatcher_threads(wf_aligner,mwavefront,score);
  const bool end_reached = wavefront_termination_end2end(wf_aligner,mwavefront,score,score_mod);
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_END_REACHED;
    wf_aligner->align_status.score = score;
    return 1; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    if (wavefront_heuristic_cufoff(wf_aligner,score,score_mod)) {
      wf_aligner->align_status.status = WF_STATUS_END_REACHED;
      return 1; // Done
    }
  }
  *max_antidiagonal = max_ak;
  return 0; // Not done
}
/*
 * Wavefront Extension (Ends-free)
 */
bool wavefront_extend_endsfree_dispatcher_seq(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int lo,
    const int hi) {
  // Parameters
  wavefront_sequences_t* const seqs = &wf_aligner->sequences;
  // Check the sequence mode
  if (seqs->mode == wf_sequences_ascii) {
    return wavefront_extend_matches_packed_endsfree(wf_aligner,mwavefront,score,lo,hi);
  } else {
    wf_offset_t dummy;
    return wavefront_extend_matches_custom(wf_aligner,mwavefront,score,lo,hi,true,&dummy);
  }
}
bool wavefront_extend_endsfree_dispatcher_threads(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score) {
  // Parameters
  const int lo = mwavefront->lo;
  const int hi = mwavefront->hi;
  bool end_reached = false;
  const int num_threads = wavefront_compute_num_threads(wf_aligner,lo,hi);
  if (num_threads == 1) {
    // Extend wavefront single-thread
    end_reached = wavefront_extend_endsfree_dispatcher_seq(wf_aligner,mwavefront,score,lo,hi);
  } else {
#ifdef WFA_PARALLEL
    // Extend wavefront in parallel
    #pragma omp parallel num_threads(num_threads)
    {
      int t_lo, t_hi;
      wavefront_compute_thread_limits(omp_get_thread_num(),omp_get_num_threads(),lo,hi,&t_lo,&t_hi);
      if (wavefront_extend_endsfree_dispatcher_seq(wf_aligner,mwavefront,score,t_lo,t_hi)) {
        end_reached = true;
      }
    }
#endif
  }
  // Return end-reached
  return end_reached;
}
int wavefront_extend_endsfree(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Modular wavefront
  const bool memory_modular = wf_aligner->wf_components.memory_modular;
  const int max_score_scope = wf_aligner->wf_components.max_score_scope;
  const int score_mod = (memory_modular) ? score % max_score_scope : score;
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score_mod];
  if (mwavefront == NULL) {
    // Check alignment feasibility (heuristic variants that can lead to no solution)
    if (wf_aligner->align_status.num_null_steps > wf_aligner->wf_components.max_score_scope) {
      wf_aligner->align_status.status = WF_STATUS_UNFEASIBLE;
      wf_aligner->align_status.score = score;
      return 1; // Done
    }
    return 0; // Not done
  }
  // Extend (dispatcher)
  const bool end_reached = wavefront_extend_endsfree_dispatcher_threads(wf_aligner,mwavefront,score);
  if (end_reached) {
    wf_aligner->align_status.status = WF_STATUS_END_REACHED;
    wf_aligner->align_status.score = score;
    return 1; // Done
  }
  // Cut-off wavefront heuristically
  if (wf_aligner->heuristic.strategy != wf_heuristic_none) {
    if (wavefront_heuristic_cufoff(wf_aligner,score,score_mod)) {
      wf_aligner->align_status.status = WF_STATUS_END_REACHED;
      return 1; // Done
    }
  }
  return 0; // Not done
}
