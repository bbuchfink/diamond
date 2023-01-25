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
 * DESCRIPTION: WaveFront alignment module for computing wavefronts (gap-affine-2p)
 */

#include "../utils/string_padded.h"
#include "wavefront_compute.h"
#include "wavefront_compute_affine.h"
#include "wavefront_backtrace_offload.h"

#ifdef WFA_PARALLEL
#include <omp.h>
#endif

/*
 * Compute Kernels
 */
void wavefront_compute_affine2p_idm(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // In Offsets
  const wf_offset_t* const m_misms = wavefront_set->in_mwavefront_misms->offsets;
  const wf_offset_t* const m_open1 = wavefront_set->in_mwavefront_open1->offsets;
  const wf_offset_t* const m_open2 = wavefront_set->in_mwavefront_open2->offsets;
  const wf_offset_t* const i1_ext = wavefront_set->in_i1wavefront_ext->offsets;
  const wf_offset_t* const i2_ext = wavefront_set->in_i2wavefront_ext->offsets;
  const wf_offset_t* const d1_ext = wavefront_set->in_d1wavefront_ext->offsets;
  const wf_offset_t* const d2_ext = wavefront_set->in_d2wavefront_ext->offsets;
  // Out Offsets
  wf_offset_t* const out_m = wavefront_set->out_mwavefront->offsets;
  wf_offset_t* const out_i1 = wavefront_set->out_i1wavefront->offsets;
  wf_offset_t* const out_i2 = wavefront_set->out_i2wavefront->offsets;
  wf_offset_t* const out_d1 = wavefront_set->out_d1wavefront->offsets;
  wf_offset_t* const out_d2 = wavefront_set->out_d2wavefront->offsets;
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE
  for (k=lo;k<=hi;++k) {
    // Update I1
    const wf_offset_t ins1_o = m_open1[k-1];
    const wf_offset_t ins1_e = i1_ext[k-1];
    const wf_offset_t ins1 = MAX(ins1_o,ins1_e) + 1;
    out_i1[k] = ins1;
    // Update I2
    const wf_offset_t ins2_o = m_open2[k-1];
    const wf_offset_t ins2_e = i2_ext[k-1];
    const wf_offset_t ins2 = MAX(ins2_o,ins2_e) + 1;
    out_i2[k] = ins2;
    // Update I
    const wf_offset_t ins = MAX(ins1,ins2);
    // Update D1
    const wf_offset_t del1_o = m_open1[k+1];
    const wf_offset_t del1_e = d1_ext[k+1];
    const wf_offset_t del1 = MAX(del1_o,del1_e);
    out_d1[k] = del1;
    // Update D2
    const wf_offset_t del2_o = m_open2[k+1];
    const wf_offset_t del2_e = d2_ext[k+1];
    const wf_offset_t del2 = MAX(del2_o,del2_e);
    out_d2[k] = del2;
    // Update D
    const wf_offset_t del = MAX(del1,del2);
    // Update M
    const wf_offset_t misms = m_misms[k] + 1;
    wf_offset_t max = MAX(del,MAX(misms,ins));
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    out_m[k] = max;
  }
}
/*
 * Compute Kernel (Piggyback)
 */
void wavefront_compute_affine2p_idm_piggyback(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  const int pattern_length = wf_aligner->pattern_length;
  const int text_length = wf_aligner->text_length;
  // In Offsets
  const wf_offset_t* const m_misms   = wavefront_set->in_mwavefront_misms->offsets;
  const wf_offset_t* const m_open1 = wavefront_set->in_mwavefront_open1->offsets;
  const wf_offset_t* const m_open2 = wavefront_set->in_mwavefront_open2->offsets;
  const wf_offset_t* const i1_ext  = wavefront_set->in_i1wavefront_ext->offsets;
  const wf_offset_t* const i2_ext  = wavefront_set->in_i2wavefront_ext->offsets;
  const wf_offset_t* const d1_ext  = wavefront_set->in_d1wavefront_ext->offsets;
  const wf_offset_t* const d2_ext  = wavefront_set->in_d2wavefront_ext->offsets;
  // Out Offsets
  wf_offset_t* const out_m  = wavefront_set->out_mwavefront->offsets;
  wf_offset_t* const out_i1 = wavefront_set->out_i1wavefront->offsets;
  wf_offset_t* const out_i2 = wavefront_set->out_i2wavefront->offsets;
  wf_offset_t* const out_d1 = wavefront_set->out_d1wavefront->offsets;
  wf_offset_t* const out_d2 = wavefront_set->out_d2wavefront->offsets;
  // In BT-pcigar
  const pcigar_t* const m_misms_bt_pcigar   = wavefront_set->in_mwavefront_misms->bt_pcigar;
  const pcigar_t* const m_open1_bt_pcigar = wavefront_set->in_mwavefront_open1->bt_pcigar;
  const pcigar_t* const m_open2_bt_pcigar = wavefront_set->in_mwavefront_open2->bt_pcigar;
  const pcigar_t* const i1_ext_bt_pcigar  = wavefront_set->in_i1wavefront_ext->bt_pcigar;
  const pcigar_t* const i2_ext_bt_pcigar  = wavefront_set->in_i2wavefront_ext->bt_pcigar;
  const pcigar_t* const d1_ext_bt_pcigar  = wavefront_set->in_d1wavefront_ext->bt_pcigar;
  const pcigar_t* const d2_ext_bt_pcigar  = wavefront_set->in_d2wavefront_ext->bt_pcigar;
  // In BT-prev
  const bt_block_idx_t* const m_misms_bt_prev   = wavefront_set->in_mwavefront_misms->bt_prev;
  const bt_block_idx_t* const m_open1_bt_prev = wavefront_set->in_mwavefront_open1->bt_prev;
  const bt_block_idx_t* const m_open2_bt_prev = wavefront_set->in_mwavefront_open2->bt_prev;
  const bt_block_idx_t* const i1_ext_bt_prev  = wavefront_set->in_i1wavefront_ext->bt_prev;
  const bt_block_idx_t* const i2_ext_bt_prev  = wavefront_set->in_i2wavefront_ext->bt_prev;
  const bt_block_idx_t* const d1_ext_bt_prev  = wavefront_set->in_d1wavefront_ext->bt_prev;
  const bt_block_idx_t* const d2_ext_bt_prev  = wavefront_set->in_d2wavefront_ext->bt_prev;
  // Out BT-pcigar
  pcigar_t* const out_m_bt_pcigar   = wavefront_set->out_mwavefront->bt_pcigar;
  pcigar_t* const out_i1_bt_pcigar  = wavefront_set->out_i1wavefront->bt_pcigar;
  pcigar_t* const out_i2_bt_pcigar  = wavefront_set->out_i2wavefront->bt_pcigar;
  pcigar_t* const out_d1_bt_pcigar  = wavefront_set->out_d1wavefront->bt_pcigar;
  pcigar_t* const out_d2_bt_pcigar  = wavefront_set->out_d2wavefront->bt_pcigar;
  // Out BT-prev
  bt_block_idx_t* const out_m_bt_prev  = wavefront_set->out_mwavefront->bt_prev;
  bt_block_idx_t* const out_i1_bt_prev = wavefront_set->out_i1wavefront->bt_prev;
  bt_block_idx_t* const out_i2_bt_prev = wavefront_set->out_i2wavefront->bt_prev;
  bt_block_idx_t* const out_d1_bt_prev = wavefront_set->out_d1wavefront->bt_prev;
  bt_block_idx_t* const out_d2_bt_prev = wavefront_set->out_d2wavefront->bt_prev;
  // Compute-Next kernel loop
  int k;
  PRAGMA_LOOP_VECTORIZE // Ifs predicated by the compiler
  for (k=lo;k<=hi;++k) {
    /*
     * I-Block
     */
    // Update I1
    const wf_offset_t ins1_o = m_open1[k-1];
    const wf_offset_t ins1_e = i1_ext[k-1];
    wf_offset_t ins1;
    pcigar_t ins1_pcigar;
    bt_block_idx_t ins1_block_idx;
    if (ins1_e >= ins1_o) {
      ins1 = ins1_e;
      ins1_pcigar = i1_ext_bt_pcigar[k-1];
      ins1_block_idx = i1_ext_bt_prev[k-1];
    } else {
      ins1 = ins1_o;
      ins1_pcigar = m_open1_bt_pcigar[k-1];
      ins1_block_idx = m_open1_bt_prev[k-1];
    }
    out_i1_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(ins1_pcigar);
    out_i1_bt_prev[k] = ins1_block_idx;
    out_i1[k] = ++ins1;
    // Update I2
    const wf_offset_t ins2_o = m_open2[k-1];
    const wf_offset_t ins2_e = i2_ext[k-1];
    wf_offset_t ins2;
    pcigar_t ins2_pcigar;
    bt_block_idx_t ins2_block_idx;
    if (ins2_e >= ins2_o) {
      ins2 = ins2_e;
      ins2_pcigar = i2_ext_bt_pcigar[k-1];
      ins2_block_idx = i2_ext_bt_prev[k-1];
    } else {
      ins2 = ins2_o;
      ins2_pcigar = m_open2_bt_pcigar[k-1];
      ins2_block_idx = m_open2_bt_prev[k-1];
    }
    out_i2_bt_pcigar[k] = PCIGAR_PUSH_BACK_INS(ins2_pcigar);
    out_i2_bt_prev[k] = ins2_block_idx;
    out_i2[k] = ++ins2;
    // Update I
    const wf_offset_t ins = MAX(ins1,ins2);
    /*
     * D-Block
     */
    // Update D1
    const wf_offset_t del1_o = m_open1[k+1];
    const wf_offset_t del1_e = d1_ext[k+1];
    wf_offset_t del1;
    pcigar_t del1_pcigar;
    bt_block_idx_t del1_block_idx;
    if (del1_e >= del1_o) {
      del1 = del1_e;
      del1_pcigar = d1_ext_bt_pcigar[k+1];
      del1_block_idx = d1_ext_bt_prev[k+1];
    } else {
      del1 = del1_o;
      del1_pcigar = m_open1_bt_pcigar[k+1];
      del1_block_idx = m_open1_bt_prev[k+1];
    }
    out_d1_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(del1_pcigar);
    out_d1_bt_prev[k] = del1_block_idx;
    out_d1[k] = del1;
    // Update D2
    const wf_offset_t del2_o = m_open2[k+1];
    const wf_offset_t del2_e = d2_ext[k+1];
    wf_offset_t del2;
    pcigar_t del2_pcigar;
    bt_block_idx_t del2_block_idx;
    if (del2_e >= del2_o) {
      del2 = del2_e;
      del2_pcigar = d2_ext_bt_pcigar[k+1];
      del2_block_idx = d2_ext_bt_prev[k+1];
    } else {
      del2 = del2_o;
      del2_pcigar = m_open2_bt_pcigar[k+1];
      del2_block_idx = m_open2_bt_prev[k+1];
    }
    out_d2_bt_pcigar[k] = PCIGAR_PUSH_BACK_DEL(del2_pcigar);
    out_d2_bt_prev[k] = del2_block_idx;
    out_d2[k] = del2;
    // Update D
    const wf_offset_t del = MAX(del1,del2);
    /*
     * M-Block
     */
    const wf_offset_t misms = m_misms[k] + 1;
    wf_offset_t max = MAX(del,MAX(misms,ins));
    if (max == ins1) {
      out_m_bt_pcigar[k] = out_i1_bt_pcigar[k];
      out_m_bt_prev[k] = out_i1_bt_prev[k];
    }
    if (max == ins2) {
      out_m_bt_pcigar[k] = out_i2_bt_pcigar[k];
      out_m_bt_prev[k] = out_i2_bt_prev[k];
    }
    if (max == del1) {
      out_m_bt_pcigar[k] = out_d1_bt_pcigar[k];
      out_m_bt_prev[k] = out_d1_bt_prev[k];
    }
    if (max == del2) {
      out_m_bt_pcigar[k] = out_d2_bt_pcigar[k];
      out_m_bt_prev[k] = out_d2_bt_prev[k];
    }
    if (max == misms) {
      out_m_bt_pcigar[k] = m_misms_bt_pcigar[k];
      out_m_bt_prev[k] = m_misms_bt_prev[k];
    }
    // Coming from I/D -> X is fake to represent gap-close
    // Coming from M -> X is real to represent mismatch
    out_m_bt_pcigar[k] = PCIGAR_PUSH_BACK_MISMS(out_m_bt_pcigar[k]);
    // Adjust offset out of boundaries !(h>tlen,v>plen) (here to allow vectorization)
    const wf_unsigned_offset_t h = WAVEFRONT_H(k,max); // Make unsigned to avoid checking negative
    const wf_unsigned_offset_t v = WAVEFRONT_V(k,max); // Make unsigned to avoid checking negative
    if (h > text_length) max = WAVEFRONT_OFFSET_NULL;
    if (v > pattern_length) max = WAVEFRONT_OFFSET_NULL;
    out_m[k] = max;
  }
}
/*
 * Compute wavefronts
 */
void wavefront_compute_affine2p_dispatcher(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  if (wavefront_set->in_mwavefront_open2->null &&
      wavefront_set->in_i2wavefront_ext->null &&
      wavefront_set->in_d2wavefront_ext->null) {
    // Delegate to regular gap-affine
    if (wf_aligner->wf_components.bt_piggyback) {
      wavefront_compute_affine_idm_piggyback(wf_aligner,wavefront_set,lo,hi);
    } else {
      wavefront_compute_affine_idm(wf_aligner,wavefront_set,lo,hi);
    }
  } else {
    // Full gap-affine-2p
    if (wf_aligner->wf_components.bt_piggyback) {
      wavefront_compute_affine2p_idm_piggyback(wf_aligner,wavefront_set,lo,hi);
    } else {
      wavefront_compute_affine2p_idm(wf_aligner,wavefront_set,lo,hi);
    }
  }
}
void wavefront_compute_affine2p_dispatcher_omp(
    wavefront_aligner_t* const wf_aligner,
    wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  const int num_threads = wavefront_compute_num_threads(wf_aligner,lo,hi);
  // Multithreading dispatcher
  if (num_threads == 1) {
    // Compute next wavefront
    wavefront_compute_affine2p_dispatcher(wf_aligner,wavefront_set,lo,hi);
  } else {
#ifdef WFA_PARALLEL
    // Compute next wavefront in parallel
    #pragma omp parallel num_threads(num_threads)
    {
      int t_lo, t_hi;
      const int thread_id = omp_get_thread_num();
      const int thread_num = omp_get_num_threads();
      wavefront_compute_thread_limits(thread_id,thread_num,lo,hi,&t_lo,&t_hi);
      wavefront_compute_affine2p_dispatcher(wf_aligner,wavefront_set,t_lo,t_hi);
    }
#endif
  }
}
void wavefront_compute_affine2p(
    wavefront_aligner_t* const wf_aligner,
    const int score) {
  // Select wavefronts
  wavefront_set_t wavefront_set;
  wavefront_compute_fetch_input(wf_aligner,&wavefront_set,score);
  // Check null wavefronts
  if (wavefront_set.in_mwavefront_misms->null &&
      wavefront_set.in_mwavefront_open1->null &&
      wavefront_set.in_mwavefront_open2->null &&
      wavefront_set.in_i1wavefront_ext->null &&
      wavefront_set.in_i2wavefront_ext->null &&
      wavefront_set.in_d1wavefront_ext->null &&
      wavefront_set.in_d2wavefront_ext->null) {
    wf_aligner->align_status.num_null_steps++; // Increment null-steps
    wavefront_compute_allocate_output_null(wf_aligner,score); // Null s-wavefront
    return;
  }
  wf_aligner->align_status.num_null_steps = 0;
  // Set limits
  int hi, lo;
  wavefront_compute_limits_input(wf_aligner,&wavefront_set,&lo,&hi);
  // Allocate wavefronts
  wavefront_compute_allocate_output(wf_aligner,&wavefront_set,score,lo,hi);
  // Init wavefront ends
  wavefront_compute_init_ends(wf_aligner,&wavefront_set,lo,hi);
  // Compute wavefronts
  wavefront_compute_affine2p_dispatcher_omp(wf_aligner,&wavefront_set,lo,hi);
  // Offload backtrace (if necessary)
  if (wf_aligner->wf_components.bt_piggyback) {
    wavefront_backtrace_offload_affine(wf_aligner,&wavefront_set,lo,hi);
  }
  // Process wavefront ends
  wavefront_compute_process_ends(wf_aligner,&wavefront_set,score);
}

