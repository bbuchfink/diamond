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
 * DESCRIPTION: WaveFront backtrace buffer to store bactrace-blocks
 */

#ifndef WAVEFRONT_BACKTRACE_BUFFER_H_
#define WAVEFRONT_BACKTRACE_BUFFER_H_

#include "../alignment/cigar.h"
#include "../utils/commons.h"
#include "../utils/vector.h"
#include "../utils/bitmap.h"
#include "../system/mm_allocator.h"
#include "wavefront_pcigar.h"
#include "wavefront_offset.h"
#include "wavefront_attributes.h"

/*
 * Separated Backtrace Block
 */
typedef uint32_t bt_block_idx_t; // Up to 2^31 references (~32GB of not-compactable pCIGARs)
#define BT_BLOCK_IDX_MAX   UINT32_MAX
#define BT_BLOCK_IDX_NULL  UINT32_MAX

typedef struct {
  pcigar_t pcigar;            // Packed CIGAR
  bt_block_idx_t prev_idx;    // Index of the previous BT-block
} __attribute__((packed)) bt_block_t;

/*
 * Backtrace initial positions
 */
typedef struct {
  int v;
  int h;
} wf_backtrace_init_pos_t;

/*
 * Backtrace Buffer
 */
typedef struct {
  // Locator
  int segment_idx;                     // Current segment idx
  int segment_offset;                  // Current free position within segment
  bt_block_t* block_next;              // Next BT-block free
  // Buffers
  vector_t* segments;                  // Memory segments (bt_block_t*)
  vector_t* alignment_init_pos;        // Buffer to store alignment's initial coordinates (h,v) (wf_backtrace_init_pos_t)
  bt_block_idx_t num_compacted_blocks; // Total compacted blocks in BT-buffer compacted (dense from 0..num_compacted_blocks-1)
  int num_compactions;                 // Total compactions performed
  // Internal buffers
  vector_t* alignment_packed;          // Temporal buffer to store final alignment (pcigar_t)
  vector_t* prefetch_blocks_idxs;      // Temporal buffer to store blocks_idxs (bt_block_idx_t)
  // MM
  mm_allocator_t* mm_allocator;
} wf_backtrace_buffer_t;

/*
 * Setup
 */
wf_backtrace_buffer_t* wf_backtrace_buffer_new(
    mm_allocator_t* const mm_allocator);
void wf_backtrace_buffer_clear(
    wf_backtrace_buffer_t* const bt_buffer);
void wf_backtrace_buffer_reap(
    wf_backtrace_buffer_t* const bt_buffer);
void wf_backtrace_buffer_delete(
    wf_backtrace_buffer_t* const bt_buffer);

/*
 * Accessors
 */
void wf_backtrace_buffer_add_used(
    wf_backtrace_buffer_t* const bt_buffer,
    const int used);
bt_block_idx_t wf_backtrace_buffer_get_mem(
    wf_backtrace_buffer_t* const bt_buffer,
    bt_block_t** const bt_block_mem,
    int* const bt_blocks_available);

/*
 * Store blocks
 */
bt_block_idx_t wf_backtrace_buffer_init_block(
    wf_backtrace_buffer_t* const bt_buffer,
    const int v,
    const int h);

/*
 * Unpack CIGAR
 */
bt_block_t* wf_backtrace_buffer_traceback_pcigar(
    wf_backtrace_buffer_t* const bt_buffer,
    bt_block_t* bt_block);
void wf_backtrace_buffer_unpack_cigar_linear(
    wf_backtrace_buffer_t* const bt_buffer,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    alignment_match_funct_t const match_funct,
    void* const match_funct_arguments,
    const int begin_v,
    const int begin_h,
    const int end_v,
    const int end_h,
    cigar_t* const cigar);
void wf_backtrace_buffer_unpack_cigar_affine(
    wf_backtrace_buffer_t* const bt_buffer,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    alignment_match_funct_t const match_funct,
    void* const match_funct_arguments,
    const int begin_v,
    const int begin_h,
    const int end_v,
    const int end_h,
    cigar_t* const cigar);

/*
 * Compact
 */
void wf_backtrace_buffer_mark_backtrace(
    wf_backtrace_buffer_t* const bt_buffer,
    const bt_block_idx_t bt_block_idx,
    bitmap_t* const bitmap);
void wf_backtrace_buffer_mark_backtrace_batch(
    wf_backtrace_buffer_t* const bt_buffer,
    wf_offset_t* const offsets,
    bt_block_idx_t* const bt_block_idxs,
    const int num_block_idxs,
    bitmap_t* const bitmap);

bt_block_idx_t wf_backtrace_buffer_compact_marked(
    wf_backtrace_buffer_t* const bt_buffer,
    bitmap_t* const bitmap,
    const int verbose);

/*
 * Utils
 */
uint64_t wf_backtrace_buffer_get_used(
    wf_backtrace_buffer_t* const bt_buffer);

bt_block_idx_t wf_backtrace_buffer_get_num_compacted_blocks(
    wf_backtrace_buffer_t* const bt_buffer);
void wf_backtrace_buffer_set_num_compacted_blocks(
    wf_backtrace_buffer_t* const bt_buffer,
    const bt_block_idx_t num_compacted_blocks);
void wf_backtrace_buffer_reset_compaction(
    wf_backtrace_buffer_t* const bt_buffer);

uint64_t wf_backtrace_buffer_get_size_allocated(
    wf_backtrace_buffer_t* const bt_buffer);
uint64_t wf_backtrace_buffer_get_size_used(
    wf_backtrace_buffer_t* const bt_buffer);

#endif /* WAVEFRONT_BACKTRACE_BUFFER_H_ */
