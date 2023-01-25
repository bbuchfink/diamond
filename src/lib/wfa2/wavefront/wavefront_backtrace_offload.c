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
 * DESCRIPTION: WaveFront alignment module for offloading partial backtraces
 */

#include "../utils/string_padded.h"
#include "wavefront_backtrace_offload.h"

/*
 * Backtrace-blocks offloading
 */
void wavefront_backtrace_offload_blocks_selective(
    wf_offset_t* const out_offsets,
    pcigar_t* const out_bt_pcigar,
    bt_block_idx_t* const out_bt_prev,
    const int lo,
    const int hi,
    const pcigar_t occupation_mask,
    wf_backtrace_buffer_t* const bt_buffer) {
  // Fetch BT-buffer free memory
  int bt_blocks_available;
  bt_block_t* bt_block_mem;
  bt_block_idx_t global_pos = wf_backtrace_buffer_get_mem(bt_buffer,&bt_block_mem,&bt_blocks_available);
  bt_block_idx_t current_pos = global_pos;
  const int max_pos = current_pos + bt_blocks_available;
  // Check PCIGAR buffers full and off-load if needed
  int k;
  for (k=lo;k<=hi;++k) {
    if (out_offsets[k]>=0 && PCIGAR_IS_UTILISED(out_bt_pcigar[k],occupation_mask)) {
      // Store
      bt_block_mem->pcigar = out_bt_pcigar[k];
      bt_block_mem->prev_idx = out_bt_prev[k];
      bt_block_mem++;
      // Reset
      out_bt_pcigar[k] = 0;
      out_bt_prev[k] = current_pos;
      current_pos++;
      // Update pos
      if (current_pos >= max_pos) {
        wf_backtrace_buffer_add_used(bt_buffer,current_pos-global_pos);
        global_pos = wf_backtrace_buffer_get_mem(bt_buffer,&bt_block_mem,&bt_blocks_available);
      }
    }
  }
  wf_backtrace_buffer_add_used(bt_buffer,current_pos-global_pos);
}
/*
 * Backtrace offloading (linear)
 */
int wavefront_backtrace_offload_blocks_linear(
    wavefront_aligner_t* const wf_aligner,
    wf_offset_t* const out_offsets,
    pcigar_t* const out_bt_pcigar,
    bt_block_idx_t* const out_bt_prev,
    const int lo,
    const int hi) {
  // Parameters
  const wavefront_memory_t wavefront_memory = wf_aligner->memory_mode;
  wf_backtrace_buffer_t* const bt_buffer = wf_aligner->wf_components.bt_buffer;
  // Select memory-mode
  switch (wavefront_memory) {
    case wavefront_memory_med:
      wavefront_backtrace_offload_blocks_selective(
          out_offsets,out_bt_pcigar,out_bt_prev,
          lo,hi,PCIGAR_HALF_FULL_MASK,bt_buffer);
      return PCIGAR_MAX_LENGTH/2; // Half occupancy
      break;
    case wavefront_memory_low:
      wavefront_backtrace_offload_blocks_selective(
          out_offsets,out_bt_pcigar,out_bt_prev,
          lo,hi,PCIGAR_FULL_MASK,bt_buffer);
      return PCIGAR_MAX_LENGTH-1; // At least 1-slots free
      break;
    default:
      fprintf(stderr,"[WFA::compute] Wrong memory-mode\n");
      exit(1);
  }
}
void wavefront_backtrace_offload_linear(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Paramters
  wavefront_t* const wf_m = wavefront_set->out_mwavefront;
  const wavefront_t* const m_misms = wavefront_set->in_mwavefront_misms;
  const wavefront_t* const m_open1 = wavefront_set->in_mwavefront_open1;
  // Compute BT occupancy maximum
  int occ_max_m = 0, occ_max_indel = 0;
  if (!m_open1->null) occ_max_indel = m_open1->bt_occupancy_max;
  if (!m_misms->null) occ_max_m = m_misms->bt_occupancy_max;
  const int occ_max = MAX(occ_max_indel,occ_max_m) + 1;
  // Set new occupancy
  wf_m->bt_occupancy_max = occ_max;
  // Offload if necessary (Gap-Linear)
  if (!wf_m->null && occ_max >= PCIGAR_MAX_LENGTH) {
    wf_offset_t* const out_m  = wavefront_set->out_mwavefront->offsets;
    pcigar_t* const out_m_bt_pcigar = wavefront_set->out_mwavefront->bt_pcigar;
    bt_block_idx_t* const out_m_bt_prev = wavefront_set->out_mwavefront->bt_prev;
    wavefront_set->out_mwavefront->bt_occupancy_max =
        wavefront_backtrace_offload_blocks_linear(
            wf_aligner,out_m,out_m_bt_pcigar,out_m_bt_prev,lo,hi);
  }
}
/*
 * Backtrace offloading (gap-affine)
 */
int wavefront_backtrace_offload_blocks_affine(
    wavefront_aligner_t* const wf_aligner,
    wf_offset_t* const out_offsets,
    pcigar_t* const out_bt_pcigar,
    bt_block_idx_t* const out_bt_prev,
    const int lo,
    const int hi) {
  // Parameters
  const wavefront_memory_t wavefront_memory = wf_aligner->memory_mode;
  wf_backtrace_buffer_t* const bt_buffer = wf_aligner->wf_components.bt_buffer;
  // Select memory-mode
  switch (wavefront_memory) {
    case wavefront_memory_med:
      wavefront_backtrace_offload_blocks_selective(
          out_offsets,out_bt_pcigar,out_bt_prev,
          lo,hi,PCIGAR_HALF_FULL_MASK,bt_buffer);
      return PCIGAR_MAX_LENGTH/2; // Half occupancy
    case wavefront_memory_low:
      wavefront_backtrace_offload_blocks_selective(
          out_offsets,out_bt_pcigar,out_bt_prev,
          lo,hi,PCIGAR_ALMOST_FULL_MASK,bt_buffer);
      return PCIGAR_MAX_LENGTH-2; // At least 2-slots free
    default:
      fprintf(stderr,"[WFA::compute] Wrong memory-mode\n");
      exit(1);
      return 0;
  }
}
void wavefront_backtrace_offload_occupation_affine(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Select distance metric
  int occ_max_m = 0;
  int occ_max_i1 = 0, occ_max_i2 = 0;
  int occ_max_d1 = 0, occ_max_d2 = 0;
  if (distance_metric == gap_affine) {
    // Parameters
    const wavefront_t* const m_misms = wavefront_set->in_mwavefront_misms;
    const wavefront_t* const m_open1 = wavefront_set->in_mwavefront_open1;
    const wavefront_t* const i1_ext = wavefront_set->in_i1wavefront_ext;
    const wavefront_t* const d1_ext = wavefront_set->in_d1wavefront_ext;
    // Compute BT occupancy maximum
    if (!m_open1->null) {
      occ_max_i1 = m_open1->bt_occupancy_max + 1;
      occ_max_d1 = m_open1->bt_occupancy_max + 1;
    }
    if (!i1_ext->null) occ_max_i1 = MAX(occ_max_i1,i1_ext->bt_occupancy_max+1);
    if (!d1_ext->null) occ_max_d1 = MAX(occ_max_d1,d1_ext->bt_occupancy_max+1);
    if (!m_misms->null) occ_max_m = m_misms->bt_occupancy_max;
    if (occ_max_m < occ_max_i1) occ_max_m = occ_max_i1;
    if (occ_max_m < occ_max_d1) occ_max_m = occ_max_d1;
    ++occ_max_m;
    // Set new occupancy
    wavefront_set->out_i1wavefront->bt_occupancy_max = occ_max_i1;
    wavefront_set->out_d1wavefront->bt_occupancy_max = occ_max_d1;
    wavefront_set->out_mwavefront->bt_occupancy_max = occ_max_m;
  } else { // distance_metric == gap_affine_2p
    // Parameters
    const wavefront_t* const m_misms = wavefront_set->in_mwavefront_misms;
    const wavefront_t* const m_open1 = wavefront_set->in_mwavefront_open1;
    const wavefront_t* const m_open2 = wavefront_set->in_mwavefront_open2;
    const wavefront_t* const i1_ext = wavefront_set->in_i1wavefront_ext;
    const wavefront_t* const i2_ext = wavefront_set->in_i2wavefront_ext;
    const wavefront_t* const d1_ext = wavefront_set->in_d1wavefront_ext;
    const wavefront_t* const d2_ext = wavefront_set->in_d2wavefront_ext;
    // Compute BT occupancy maximum (I)
    if (!m_open1->null) {
      occ_max_i1 = m_open1->bt_occupancy_max + 1;
      occ_max_d1 = m_open1->bt_occupancy_max + 1;
    }
    // Compute BT occupancy maximum (D)
    if (!i1_ext->null) occ_max_i1 = MAX(occ_max_i1,i1_ext->bt_occupancy_max+1);
    if (!d1_ext->null) occ_max_d1 = MAX(occ_max_d1,d1_ext->bt_occupancy_max+1);
    if (!m_open2->null) {
      occ_max_i2 = m_open2->bt_occupancy_max + 1;
      occ_max_d2 = m_open2->bt_occupancy_max + 1;
    }
    if (!i2_ext->null) occ_max_i2 = MAX(occ_max_i2,i2_ext->bt_occupancy_max+1);
    if (!d2_ext->null) occ_max_d2 = MAX(occ_max_d2,d2_ext->bt_occupancy_max+1);
    // Compute BT occupancy maximum (M)
    if (!m_misms->null) occ_max_m = m_misms->bt_occupancy_max;
    if (occ_max_m < occ_max_i1) occ_max_m = occ_max_i1;
    if (occ_max_m < occ_max_i2) occ_max_m = occ_max_i2;
    if (occ_max_m < occ_max_d1) occ_max_m = occ_max_d1;
    if (occ_max_m < occ_max_d2) occ_max_m = occ_max_d2;
    ++occ_max_m;
    // Set new occupancy
    wavefront_set->out_i1wavefront->bt_occupancy_max = occ_max_i1;
    wavefront_set->out_i2wavefront->bt_occupancy_max = occ_max_i2;
    wavefront_set->out_d1wavefront->bt_occupancy_max = occ_max_d1;
    wavefront_set->out_d2wavefront->bt_occupancy_max = occ_max_d2;
    wavefront_set->out_mwavefront->bt_occupancy_max = occ_max_m;
  }
}
void wavefront_backtrace_offload_affine(
    wavefront_aligner_t* const wf_aligner,
    const wavefront_set_t* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->penalties.distance_metric;
  // Compute maximum occupancy
  wavefront_backtrace_offload_occupation_affine(wf_aligner,wavefront_set);
  // Offload if necessary (Gap-Affine)
  const wavefront_t* const wf_m = wavefront_set->out_mwavefront;
  if (!wf_m->null && wf_m->bt_occupancy_max >= PCIGAR_MAX_LENGTH-1) {
    wf_offset_t* const out_m  = wavefront_set->out_mwavefront->offsets;
    pcigar_t* const out_m_bt_pcigar = wavefront_set->out_mwavefront->bt_pcigar;
    bt_block_idx_t* const out_m_bt_prev = wavefront_set->out_mwavefront->bt_prev;
    wavefront_set->out_mwavefront->bt_occupancy_max =
        wavefront_backtrace_offload_blocks_affine(
            wf_aligner,out_m,out_m_bt_pcigar,out_m_bt_prev,lo,hi);
  }
  const wavefront_t* const wf_i1 = wavefront_set->out_i1wavefront;
  if (!wf_i1->null && wf_i1->bt_occupancy_max >= PCIGAR_MAX_LENGTH-1) {
    wf_offset_t* const out_i1 = wavefront_set->out_i1wavefront->offsets;
    pcigar_t* const out_i1_bt_pcigar = wavefront_set->out_i1wavefront->bt_pcigar;
    bt_block_idx_t* const out_i1_bt_prev = wavefront_set->out_i1wavefront->bt_prev;
    wavefront_set->out_i1wavefront->bt_occupancy_max =
        wavefront_backtrace_offload_blocks_affine(
            wf_aligner,out_i1,out_i1_bt_pcigar,out_i1_bt_prev,lo,hi);
  }
  const wavefront_t* const wf_d1 = wavefront_set->out_d1wavefront;
  if (!wf_d1->null && wf_d1->bt_occupancy_max >= PCIGAR_MAX_LENGTH-1) {
    wf_offset_t* const out_d1 = wavefront_set->out_d1wavefront->offsets;
    pcigar_t* const out_d1_bt_pcigar  = wavefront_set->out_d1wavefront->bt_pcigar;
    bt_block_idx_t* const out_d1_bt_prev = wavefront_set->out_d1wavefront->bt_prev;
    wavefront_set->out_d1wavefront->bt_occupancy_max =
        wavefront_backtrace_offload_blocks_affine(
            wf_aligner,out_d1,out_d1_bt_pcigar,out_d1_bt_prev,lo,hi);
  }
  if (distance_metric == gap_affine) return;
  // Offload if necessary (Gap-Affine-2p)
  const wavefront_t* const wf_i2 = wavefront_set->out_i2wavefront;
  if (!wf_i2->null && wf_i2->bt_occupancy_max >= PCIGAR_MAX_LENGTH-1) {
    wf_offset_t* const out_i2 = wavefront_set->out_i2wavefront->offsets;
    pcigar_t* const out_i2_bt_pcigar = wavefront_set->out_i2wavefront->bt_pcigar;
    bt_block_idx_t* const out_i2_bt_prev = wavefront_set->out_i2wavefront->bt_prev;
    wavefront_set->out_i2wavefront->bt_occupancy_max =
        wavefront_backtrace_offload_blocks_affine(
            wf_aligner,out_i2,out_i2_bt_pcigar,out_i2_bt_prev,lo,hi);
  }
  const wavefront_t* const wf_d2 = wavefront_set->out_d2wavefront;
  if (!wf_d2->null && wf_d2->bt_occupancy_max >= PCIGAR_MAX_LENGTH-1) {
    wf_offset_t* const out_d2 = wavefront_set->out_d2wavefront->offsets;
    pcigar_t* const out_d2_bt_pcigar = wavefront_set->out_d2wavefront->bt_pcigar;
    bt_block_idx_t* const out_d2_bt_prev = wavefront_set->out_d2wavefront->bt_prev;
    wavefront_set->out_d2wavefront->bt_occupancy_max =
        wavefront_backtrace_offload_blocks_affine(
            wf_aligner,out_d2,out_d2_bt_pcigar,out_d2_bt_prev,lo,hi);
  }
}

