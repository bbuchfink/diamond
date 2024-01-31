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

#include "utils/commons.h"
#include "system/mm_allocator.h"
#include "wavefront_unialign.h"
#include "wavefront.h"
#include "wavefront_attributes.h"
#include "wavefront_offset.h"
#include "wavefront_penalties.h"
#include "wavefront_plot.h"
#include "wavefront_slab.h"

#include "wavefront_components.h"
#include "wavefront_compute.h"
#include "wavefront_compute_affine.h"
#include "wavefront_compute_affine2p.h"
#include "wavefront_compute_edit.h"
#include "wavefront_compute_linear.h"
#include "wavefront_extend.h"
#include "wavefront_backtrace.h"
#include "wavefront_backtrace_buffer.h"

/*
 * Initialize alignment
 */
void wavefront_unialign_init(
    wavefront_aligner_t* const wf_aligner,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end) {
  // Parameters
  wavefront_align_status_t* const align_status = &wf_aligner->align_status;
  alignment_form_t* const alignment_form = &wf_aligner->alignment_form;
  const bool end2end = (alignment_form->span == alignment_end2end);
  // Configure WF-compute function
  switch (wf_aligner->penalties.distance_metric) {
    case indel:
    case edit:
      align_status->wf_align_compute = &wavefront_compute_edit;
      break;
    case gap_linear:
      align_status->wf_align_compute = &wavefront_compute_linear;
      break;
    case gap_affine:
      align_status->wf_align_compute = &wavefront_compute_affine;
      break;
    case gap_affine_2p:
      align_status->wf_align_compute = &wavefront_compute_affine2p;
      break;
    default:
      fprintf(stderr,"[WFA] Distance function not implemented\n");
      exit(1);
      break;
  }
  // Configure WF-extend function
  if (end2end) {
    align_status->wf_align_extend = &wavefront_extend_end2end;
  } else {
    align_status->wf_align_extend = &wavefront_extend_endsfree;
  }
  // Initialize wavefront-aligner (to perform a new alignment)
  wf_aligner->component_begin = component_begin;
  wf_aligner->component_end = component_end;
  wavefront_aligner_init(wf_aligner,0);
  // Clear cigar
  cigar_clear(wf_aligner->cigar);
}
/*
 * Limits
 */
bool wavefront_unialign_reached_limits(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Check alignment-score limit
  if (score >= wf_aligner->system.max_alignment_steps) {
    wf_aligner->cigar->score = -wf_aligner->system.max_alignment_steps;
    wf_aligner->align_status.status = WF_STATUS_MAX_STEPS_REACHED;
    wf_aligner->align_status.score = score;
    return true; // Stop
  }
  // Global probing interval
  alignment_system_t* const system = &wf_aligner->system;
  if (score % system->probe_interval_global != 0) return false; // Continue
  if (system->verbose >= 3) {
    wavefront_unialign_print_status(stderr,wf_aligner,score); // DEBUG
  }
  // BT-Buffer
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  if (wf_components->bt_buffer!=NULL && (score%system->probe_interval_compact)==0) {
    uint64_t bt_memory = wf_backtrace_buffer_get_size_used(wf_components->bt_buffer);
    // Check BT-buffer memory
    if (bt_memory > system->max_memory_compact) {
      // Compact BT-buffer
      wavefront_components_compact_bt_buffer(wf_components,score,wf_aligner->system.verbose);
      // Set new buffer limit
      bt_memory = wf_backtrace_buffer_get_size_used(wf_components->bt_buffer);
      uint64_t proposed_mem = (double)bt_memory * TELESCOPIC_FACTOR;
      if (system->max_memory_compact < proposed_mem && proposed_mem < system->max_memory_abort) {
        proposed_mem = system->max_memory_compact;
      }
      // Reset (if maximum compacts has been performed)
      if (wf_components->bt_buffer->num_compactions >= system->max_partial_compacts) {
        wf_backtrace_buffer_reset_compaction(wf_components->bt_buffer);
      }
    }
  }
  // Check overall memory used
  const uint64_t wf_memory_used = wavefront_aligner_get_size(wf_aligner);
  if (wf_memory_used > system->max_memory_abort) {
    wf_aligner->align_status.status = WF_STATUS_OOM;
    wf_aligner->align_status.score = score;
    return true; // Stop
  }
  // Otherwise continue
  return false;
}
/*
 * Terminate alignment (backtrace)
 */
void wavefront_unialign_terminate(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Parameters
  wavefront_align_status_t* const align_status = &wf_aligner->align_status;
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  cigar_t* const cigar = wf_aligner->cigar;
  // Select alignment scope
  align_status->score = score;
  if (wf_aligner->alignment_scope == compute_score) {
    // Set end-alignment position & score
    if (align_status->status == WF_STATUS_END_REACHED) {
      cigar->end_v = pattern_length;
      cigar->end_h = text_length;
      cigar->score = wavefront_compute_classic_score(wf_aligner,pattern_length,text_length,score);
      align_status->status = WF_STATUS_ALG_COMPLETED;
    } else {
      const int k = wf_aligner->alignment_end_pos.k;
      const int offset = wf_aligner->alignment_end_pos.offset;
      cigar->end_v = WAVEFRONT_V(k,offset);
      cigar->end_h = WAVEFRONT_H(k,offset);
      cigar->score = wavefront_compute_classic_score(wf_aligner,cigar->end_v,cigar->end_h,score);
      align_status->dropped = true;
      align_status->status = WF_STATUS_ALG_PARTIAL;
    }
  } else {
    // Parameters
    wavefront_components_t* const wf_components = &wf_aligner->wf_components;
    const int alignment_end_k = wf_aligner->alignment_end_pos.k;
    const wf_offset_t alignment_end_offset = wf_aligner->alignment_end_pos.offset;
    if (alignment_end_offset != WAVEFRONT_OFFSET_NULL) {
      if (wf_components->bt_piggyback) {
        // Fetch wavefront
        const bool memory_modular = wf_aligner->wf_components.memory_modular;
        const int max_score_scope = wf_aligner->wf_components.max_score_scope;
        const int score_mod = (memory_modular) ? score % max_score_scope : score;
        wavefront_t* const mwavefront = wf_components->mwavefronts[score_mod];
        // Backtrace alignment from buffer (unpacking pcigar)
        wavefront_backtrace_pcigar(
            wf_aligner,alignment_end_k,alignment_end_offset,
            mwavefront->bt_pcigar[alignment_end_k],
            mwavefront->bt_prev[alignment_end_k]);
      } else {
        // Backtrace alignment
        if (wf_aligner->penalties.distance_metric <= gap_linear) {
          wavefront_backtrace_linear(wf_aligner,
              score,alignment_end_k,alignment_end_offset);
        } else {
          wavefront_backtrace_affine(wf_aligner,
              wf_aligner->component_begin,wf_aligner->component_end,
              score,alignment_end_k,alignment_end_offset);
        }
      }
    }
    /*
     * Post-processing (Extension-Trim, Score, and Ends)
     *
     *                   |     Alignment-Regular    |      Alignment-Extension         |
     *  |------------------------------------------------------------------------------|
     *  |  END_REACHABLE |  NoTrim + ALG_COMPLETED  | Trim + ALG_PARTIAL|ALG_COMPLETED |
     *  |END_UNREACHABLE |  Trim + ALG_PARTIAL      | Trim + ALG_PARTIAL               |
     */
    const bool do_extension = wf_aligner->alignment_form.extension;
    const bool unreachable = (align_status->status == WF_STATUS_END_UNREACHABLE);
    align_status->dropped = unreachable;
    if (do_extension || unreachable) {
      // Alignment extension (maximal score)
      const bool cigar_trimmed = wavefront_aligner_maxtrim_cigar(wf_aligner);
      if (cigar_trimmed) {
        align_status->status = WF_STATUS_ALG_PARTIAL;
      } else {
        align_status->status = (align_status->status == WF_STATUS_END_UNREACHABLE) ?
            WF_STATUS_ALG_PARTIAL : WF_STATUS_ALG_COMPLETED;
      }
    } else {
      const int k = wf_aligner->alignment_end_pos.k;
      const int offset = wf_aligner->alignment_end_pos.offset;
      cigar->end_v = WAVEFRONT_V(k,offset);
      cigar->end_h = WAVEFRONT_H(k,offset);
      cigar->score = wavefront_compute_classic_score(wf_aligner,cigar->end_v,cigar->end_h,score);
      // Set status
      if (unreachable) {
        align_status->status = WF_STATUS_ALG_PARTIAL;
      } else {
        align_status->status = WF_STATUS_ALG_COMPLETED;
      }
    }
  }
}
/*
 * Classic WF-Alignment (Unidirectional)
 */
#define WFA_UNIALIGN_DEBUG() wavefront_aligner_print(stderr,wf_aligner,0,score,7,0)
int wavefront_unialign(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  wavefront_align_status_t* const align_status = &wf_aligner->align_status;
  void (*wf_align_compute)(wavefront_aligner_t* const,const int) = align_status->wf_align_compute;
  int (*wf_align_extend)(wavefront_aligner_t* const,const int) = align_status->wf_align_extend;
  // Compute wavefronts of increasing score
  int score = align_status->score;
  // WFA_UNIALIGN_DEBUG(); // DEBUG
  while (true) {
    // Exact extend s-wavefront
    const int finished = (*wf_align_extend)(wf_aligner,score);
    if (finished) {
      // WFA_UNIALIGN_DEBUG(); // DEBUG
      if (align_status->status == WF_STATUS_END_REACHED ||
          align_status->status == WF_STATUS_END_UNREACHABLE) {
        wavefront_unialign_terminate(wf_aligner,score);
      }
      return align_status->status;
    }
    // Compute (s+1)-wavefront
    ++score;
    (*wf_align_compute)(wf_aligner,score);
    // Probe limits
    if (wavefront_unialign_reached_limits(wf_aligner,score)) return align_status->status;
    // Plot
    if (wf_aligner->plot != NULL) wavefront_plot(wf_aligner,score,0);
    // WFA_UNIALIGN_DEBUG(); // DEBUG
  }
  // Unreachable code
  return WF_STATUS_OK;
}
/*
 * Display
 */
void wavefront_unialign_print_status(
    FILE* const stream,
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Parameters
  wavefront_components_t* const wf_components = &wf_aligner->wf_components;
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  // Approximate progress
  const int dist_total = MAX(text_length,pattern_length);
  int s = (wf_components->memory_modular) ? score%wf_components->max_score_scope : score;
  wavefront_t* wavefront = wf_components->mwavefronts[s];
  if (wavefront==NULL && s>0) {
    s = (wf_components->memory_modular) ? (score-1)%wf_components->max_score_scope : (score-1);
    wavefront = wf_components->mwavefronts[s];
  }
  int dist_max = -1, wf_len = -1, k;
  if (wavefront!=NULL) {
    wf_offset_t* const offsets = wavefront->offsets;
    for (k=wavefront->lo;k<=wavefront->hi;++k) {
      const int dist = MAX(WAVEFRONT_V(k,offsets[k]),WAVEFRONT_H(k,offsets[k]));
      dist_max = MAX(dist_max,dist);
    }
    wf_len = wavefront->hi-wavefront->lo+1;
  }
  // Memory used
  const uint64_t slab_size = wavefront_slab_get_size(wf_aligner->wavefront_slab);
  const uint64_t bt_buffer_used = (wf_components->bt_buffer) ?
      wf_backtrace_buffer_get_size_used(wf_components->bt_buffer) : 0;
  // Progress
  const float aligned_progress = (dist_max>=0) ? (100.0f*(float)dist_max/(float)dist_total) : -1.0f;
  const float million_offsets = (wf_len>=0) ? (float)wf_len/1000000.0f : -1.0f;
  // Print one-line status
  fprintf(stream,"[");
  wavefront_aligner_print_mode(stream,wf_aligner);
  fprintf(stream,
      "] SequenceLength=(%d,%d) Score %d (~ %2.3f%% aligned). "
      "MemoryUsed(WF-Slab,BT-buffer)=(%lu MB,%lu MB). "
      "Wavefronts ~ %2.3f Moffsets\n",
      pattern_length,
      text_length,
      score,
      aligned_progress,
      CONVERT_B_TO_MB(slab_size),
      CONVERT_B_TO_MB(bt_buffer_used),
      million_offsets);
}
