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

#include "wavefront_backtrace_buffer.h"

/*
 * Config
 */
#define BT_BUFFER_SEGMENT_LENGTH BUFFER_SIZE_8M

#define BT_BUFFER_SEGMENT_IDX(block_idx)    ((block_idx)/BT_BUFFER_SEGMENT_LENGTH)
#define BT_BUFFER_SEGMENT_OFFSET(block_idx) ((block_idx)%BT_BUFFER_SEGMENT_LENGTH)

#define BT_BUFFER_IDX(segment_idx,segment_offset) \
  ((segment_idx)*BT_BUFFER_SEGMENT_LENGTH) + (segment_offset)

/*
 * BT-Block Segments
 */
void wf_backtrace_buffer_segment_add(
    wf_backtrace_buffer_t* const bt_buffer) {
  bt_block_t* const bt_segment = mm_allocator_calloc(
      bt_buffer->mm_allocator,BT_BUFFER_SEGMENT_LENGTH,bt_block_t,false);
  vector_insert(bt_buffer->segments,bt_segment,bt_block_t*);
}
void wf_backtrace_buffer_segment_reserve(
    wf_backtrace_buffer_t* const bt_buffer) {
  // Reset position
  bt_buffer->segment_offset = 0;
  ++(bt_buffer->segment_idx);
  // Check segments
  if (bt_buffer->segment_idx >= vector_get_used(bt_buffer->segments)) {
    // Check segment position
    const uint64_t block_idx = ((uint64_t)bt_buffer->segment_idx+1) * BT_BUFFER_SEGMENT_LENGTH;
    if (block_idx >= BT_BLOCK_IDX_MAX) {
      fprintf(stderr,"[WFA::BacktraceBuffer] Reached maximum addressable index"); exit(-1);
    }
    // Add segment
    wf_backtrace_buffer_segment_add(bt_buffer);
  }
  // Set pointer to next block free
  bt_block_t** const segments = vector_get_mem(bt_buffer->segments,bt_block_t*);
  bt_buffer->block_next = segments[bt_buffer->segment_idx];
}
/*
 * Setup
 */
wf_backtrace_buffer_t* wf_backtrace_buffer_new(
    mm_allocator_t* const mm_allocator) {
  // Alloc
  wf_backtrace_buffer_t* const bt_buffer =
      mm_allocator_alloc(mm_allocator,wf_backtrace_buffer_t);
  bt_buffer->mm_allocator = mm_allocator;
  // Initialize
  bt_buffer->segment_idx = 0;
  bt_buffer->segment_offset = 0;
  bt_buffer->segments = vector_new(10,bt_block_t*);
  wf_backtrace_buffer_segment_add(bt_buffer); // Add initial segment
  bt_buffer->block_next = vector_get_mem(bt_buffer->segments,bt_block_t*)[0];
  bt_buffer->num_compacted_blocks = 0;
  bt_buffer->num_compactions = 0;
  bt_buffer->alignment_init_pos = vector_new(100,wf_backtrace_init_pos_t);
  bt_buffer->alignment_packed = vector_new(100,pcigar_t);
  bt_buffer->prefetch_blocks_idxs = vector_new(500,bt_block_idx_t);
  // Return
  return bt_buffer;
}
void wf_backtrace_buffer_clear(
    wf_backtrace_buffer_t* const bt_buffer) {
  bt_buffer->segment_idx = 0;
  bt_buffer->segment_offset = 0;
  bt_buffer->block_next = vector_get_mem(bt_buffer->segments,bt_block_t*)[0];
  bt_buffer->num_compacted_blocks = 0;
  bt_buffer->num_compactions = 0;
  vector_clear(bt_buffer->alignment_init_pos);
}
void wf_backtrace_buffer_reap(
    wf_backtrace_buffer_t* const bt_buffer) {
  // Reap segments beyond the first
  const int num_segments = vector_get_used(bt_buffer->segments);
  bt_block_t** const segments = vector_get_mem(bt_buffer->segments,bt_block_t*);
  int i;
  for (i=1;i<num_segments;++i) {
    mm_allocator_free(bt_buffer->mm_allocator,segments[i]);
  }
  vector_set_used(bt_buffer->segments,1);
  // Clear
  bt_buffer->segment_idx = 0;
  bt_buffer->segment_offset = 0;
  bt_buffer->block_next = vector_get_mem(bt_buffer->segments,bt_block_t*)[0];
  bt_buffer->num_compacted_blocks = 0;
  bt_buffer->num_compactions = 0;
}
void wf_backtrace_buffer_delete(
    wf_backtrace_buffer_t* const bt_buffer) {
  // Free segments
  const int num_segments = vector_get_used(bt_buffer->segments);
  bt_block_t** const segments = vector_get_mem(bt_buffer->segments,bt_block_t*);
  int i;
  for (i=0;i<num_segments;++i) {
    mm_allocator_free(bt_buffer->mm_allocator,segments[i]);
  }
  // Free handlers
  vector_delete(bt_buffer->segments);
  vector_delete(bt_buffer->alignment_init_pos);
  vector_delete(bt_buffer->alignment_packed);
  vector_delete(bt_buffer->prefetch_blocks_idxs);
  mm_allocator_free(bt_buffer->mm_allocator,bt_buffer);
}
/*
 * Accessors
 */
uint64_t wf_backtrace_buffer_get_used(
    wf_backtrace_buffer_t* const bt_buffer) {
  const bt_block_idx_t max_block_idx = BT_BUFFER_IDX(bt_buffer->segment_idx,bt_buffer->segment_offset);
  return max_block_idx;
}
bt_block_idx_t wf_backtrace_buffer_get_num_compacted_blocks(
    wf_backtrace_buffer_t* const bt_buffer) {
  return bt_buffer->num_compacted_blocks;
}
void wf_backtrace_buffer_set_num_compacted_blocks(
    wf_backtrace_buffer_t* const bt_buffer,
    const bt_block_idx_t num_compacted_blocks) {
  bt_buffer->num_compacted_blocks = num_compacted_blocks;
}
void wf_backtrace_buffer_reset_compaction(
    wf_backtrace_buffer_t* const bt_buffer) {
  bt_buffer->num_compactions = 0;
  bt_buffer->num_compacted_blocks = 0;
}
uint64_t wf_backtrace_buffer_get_size_allocated(
    wf_backtrace_buffer_t* const bt_buffer) {
  const uint64_t segments_used = vector_get_used(bt_buffer->segments);
  return segments_used*BT_BUFFER_SEGMENT_LENGTH*sizeof(bt_block_t);
}
uint64_t wf_backtrace_buffer_get_size_used(
    wf_backtrace_buffer_t* const bt_buffer) {
  const bt_block_idx_t max_block_idx = BT_BUFFER_IDX(bt_buffer->segment_idx,bt_buffer->segment_offset);
  return max_block_idx*sizeof(bt_block_t);
}
void wf_backtrace_buffer_prefetch_block(
    wf_backtrace_buffer_t* const bt_buffer,
    const bt_block_idx_t block_idx) {
  // Compute location
  const int segment_idx = BT_BUFFER_SEGMENT_IDX(block_idx);
  const int segment_offset = BT_BUFFER_SEGMENT_OFFSET(block_idx);
  // Fetch bt-block
  bt_block_t** const segments = vector_get_mem(bt_buffer->segments,bt_block_t*);
  PREFETCH(segments[segment_idx]+segment_offset);
}
bt_block_t* wf_backtrace_buffer_get_block(
    wf_backtrace_buffer_t* const bt_buffer,
    const bt_block_idx_t block_idx) {
  // Compute location
  const int segment_idx = BT_BUFFER_SEGMENT_IDX(block_idx);
  const int segment_offset = BT_BUFFER_SEGMENT_OFFSET(block_idx);
  // Fetch bt-block
  bt_block_t** const segments = vector_get_mem(bt_buffer->segments,bt_block_t*);
  return &(segments[segment_idx][segment_offset]);
}
void wf_backtrace_buffer_add_used(
    wf_backtrace_buffer_t* const bt_buffer,
    const int used) {
  // Next
  bt_buffer->segment_offset += used;
  bt_buffer->block_next += used;
  // Reserve
  if (bt_buffer->segment_offset >= BT_BUFFER_SEGMENT_LENGTH) {
    wf_backtrace_buffer_segment_reserve(bt_buffer);
  }
}
bt_block_idx_t wf_backtrace_buffer_get_mem(
    wf_backtrace_buffer_t* const bt_buffer,
    bt_block_t** const bt_block_mem,
    int* const bt_blocks_available) {
  // Parameters
  const int segment_idx = bt_buffer->segment_idx;
  const int segment_offset = bt_buffer->segment_offset;
  // Get total available blocks
  *bt_block_mem = bt_buffer->block_next;
  *bt_blocks_available = BT_BUFFER_SEGMENT_LENGTH - bt_buffer->segment_offset;
  // Return current global position
  return BT_BUFFER_IDX(segment_idx,segment_offset);
}
/*
 * Store blocks
 */
void wf_backtrace_buffer_store_block(
    wf_backtrace_buffer_t* const bt_buffer,
    const pcigar_t pcigar,
    const bt_block_idx_t prev_idx) {
  // Store BT-block
  bt_buffer->block_next->pcigar = pcigar;
  bt_buffer->block_next->prev_idx = prev_idx;
  // Next
  ++(bt_buffer->block_next);
  ++(bt_buffer->segment_offset);
  // Reserve
  if (bt_buffer->segment_offset >= BT_BUFFER_SEGMENT_LENGTH) {
    wf_backtrace_buffer_segment_reserve(bt_buffer);
  }
}
bt_block_idx_t wf_backtrace_buffer_init_block(
    wf_backtrace_buffer_t* const bt_buffer,
    const int v,
    const int h) {
  // Parameters
  const int segment_idx = bt_buffer->segment_idx;
  const int segment_offset = bt_buffer->segment_offset;
  // Store initial position (v,h)
  const int init_position_offset = vector_get_used(bt_buffer->alignment_init_pos);
  wf_backtrace_init_pos_t init_pos = { .v = v, .h = h };
  vector_insert(bt_buffer->alignment_init_pos,init_pos,wf_backtrace_init_pos_t);
  // Store BT-block (Index to initial position,NULL prev)
  wf_backtrace_buffer_store_block(bt_buffer,init_position_offset,BT_BLOCK_IDX_NULL);
  // Return current index
  return BT_BUFFER_IDX(segment_idx,segment_offset);
}
/*
 * Unpack CIGAR
 */
bt_block_t* wf_backtrace_buffer_traceback_pcigar(
    wf_backtrace_buffer_t* const bt_buffer,
    bt_block_t* bt_block) {
  // Clear temporal buffer
  vector_t* const alignment_packed = bt_buffer->alignment_packed;
  vector_clear(alignment_packed);
  // Traverse-back the BT-blocks and store all the pcigars
  while (bt_block->prev_idx != BT_BLOCK_IDX_NULL) {
    vector_insert(alignment_packed,bt_block->pcigar,pcigar_t);
    const bt_block_idx_t prev_idx = bt_block->prev_idx;
    bt_block = wf_backtrace_buffer_get_block(bt_buffer,prev_idx);
  }
  // Return initial block (start coordinate)
  return bt_block;
}
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
    cigar_t* const cigar) {
  // Clear cigar
  char* cigar_buffer = cigar->operations;
  cigar->begin_offset = 0;
  // Add init insertions/deletions
  int i;
  int v = begin_v;
  int h = begin_h;
  for (i=0;i<h;++i) {*cigar_buffer = 'I'; ++cigar_buffer;};
  for (i=0;i<v;++i) {*cigar_buffer = 'D'; ++cigar_buffer;};
  // Traverse-forward the pcigars and unpack the cigar
  const int num_palignment_blocks = vector_get_used(bt_buffer->alignment_packed);
  pcigar_t* const palignment_blocks = vector_get_mem(bt_buffer->alignment_packed,pcigar_t);
  for (i=num_palignment_blocks-1;i>=0;--i) {
    // Unpack block
    int cigar_block_length = 0;
    pcigar_unpack_linear(
        palignment_blocks[i],
        pattern,pattern_length,text,text_length,
        match_funct,match_funct_arguments,&v,&h,
        cigar_buffer,&cigar_block_length);
    // Update CIGAR
    cigar_buffer += cigar_block_length;
  }
  // Account for last stroke of matches
  const int num_matches = MIN(end_v-v,end_h-h);
  for (i=0;i<num_matches;++i) {*cigar_buffer = 'M'; ++cigar_buffer;};
  v += num_matches;
  h += num_matches;
  // Account for last stroke of insertion/deletion
  while (h < text_length) {*cigar_buffer = 'I'; ++cigar_buffer; ++h;};
  while (v < pattern_length) {*cigar_buffer = 'D'; ++cigar_buffer; ++v;};
  // Close CIGAR
  *cigar_buffer = '\0';
  cigar->end_offset = cigar_buffer - cigar->operations;
}
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
    cigar_t* const cigar) {
  // Clear cigar
  char* cigar_buffer = cigar->operations;
  cigar->begin_offset = 0;
  // Add init insertions/deletions
  int i;
  int v = begin_v;
  int h = begin_h;
  for (i=0;i<h;++i) {*cigar_buffer = 'I'; ++cigar_buffer;};
  for (i=0;i<v;++i) {*cigar_buffer = 'D'; ++cigar_buffer;};
  // Traverse-forward the pcigars and unpack the cigar
  const int num_palignment_blocks = vector_get_used(bt_buffer->alignment_packed);
  pcigar_t* const palignment_blocks = vector_get_mem(bt_buffer->alignment_packed,pcigar_t);
  affine_matrix_type current_matrix_type = affine_matrix_M;
  for (i=num_palignment_blocks-1;i>=0;--i) {
    // Unpack block
    int cigar_block_length = 0;
    pcigar_unpack_affine(
        palignment_blocks[i],
        pattern,pattern_length,text,text_length,
        match_funct,match_funct_arguments,&v,&h,
        cigar_buffer,&cigar_block_length,&current_matrix_type);
    // Update CIGAR
    cigar_buffer += cigar_block_length;
  }
  // Account for last stroke of matches
  const int num_matches = MIN(end_v-v,end_h-h);
  for (i=0;i<num_matches;++i) {*cigar_buffer = 'M'; ++cigar_buffer;};
  v += num_matches;
  h += num_matches;
  // Account for last stroke of insertion/deletion
  while (h < text_length) {*cigar_buffer = 'I'; ++cigar_buffer; ++h;};
  while (v < pattern_length) {*cigar_buffer = 'D'; ++cigar_buffer; ++v;};
  // Close CIGAR
  *cigar_buffer = '\0';
  cigar->end_offset = cigar_buffer - cigar->operations;
}
/*
 * Compact
 */
void wf_backtrace_buffer_mark_backtrace(
    wf_backtrace_buffer_t* const bt_buffer,
    const bt_block_idx_t bt_block_idx,
    bitmap_t* const bitmap) {
  // Parameters
  const bt_block_idx_t num_compacted_blocks = bt_buffer->num_compacted_blocks;
  // Traverse-back the BT-blocks while not marked
  bt_block_t bt_block_last = { .prev_idx = bt_block_idx };
  bt_block_t* bt_block = &bt_block_last;
  // Check marked and fetch previous (until already marked or NULL is found)
  while (bt_block->prev_idx != BT_BLOCK_IDX_NULL &&
         bt_block->prev_idx >= num_compacted_blocks &&
         !bitmap_check__set(bitmap,bt_block->prev_idx)) {
    // Fetch previous BT-block
    const bt_block_idx_t prev_idx = bt_block->prev_idx;
    bt_block = wf_backtrace_buffer_get_block(bt_buffer,prev_idx);
  }
}
void wf_backtrace_buffer_mark_backtrace_batch(
    wf_backtrace_buffer_t* const bt_buffer,
    wf_offset_t* const offsets,
    bt_block_idx_t* const bt_block_idxs,
    const int num_block_idxs,
    bitmap_t* const bitmap) {
  // Parameters
  const bt_block_idx_t num_compacted_blocks = bt_buffer->num_compacted_blocks;
  // Reserve prefetch-buffer
  const int max_batch_size = 100;
  vector_reserve(bt_buffer->prefetch_blocks_idxs,max_batch_size,false);
  bt_block_idx_t* const pf_block_idx = vector_get_mem(bt_buffer->prefetch_blocks_idxs,bt_block_idx_t);
  // Fill-in loop (+ initial prefetch)
  int active_blocks = 0, next_idx = 0;
  while (active_blocks < max_batch_size && next_idx < num_block_idxs) {
    // Check NULL
    const bt_block_idx_t block_idx = bt_block_idxs[next_idx];
    if (offsets[next_idx] >= 0 && 
        block_idx >= num_compacted_blocks) { // NOTE block_idx != BT_BLOCK_IDX_NULL
      // Prefetch (bt-block and bt_block)
      BITMAP_PREFETCH_BLOCK(bitmap,block_idx);
      wf_backtrace_buffer_prefetch_block(bt_buffer,block_idx);
      // Store
      pf_block_idx[active_blocks] = block_idx;
      ++active_blocks;
    }
    ++next_idx; // Next
  }
  // Batch process+prefetch loop
  int i = 0;
  while (active_blocks > 0) {
    // Fetch BT-block & BM-block
    const bt_block_idx_t block_idx = pf_block_idx[i];
    BITMAP_GET_BLOCK(bitmap,block_idx,block_bm_ptr);
    // Check marked
    if (!BM_BLOCK_IS_SET(*block_bm_ptr,block_idx)) {
      BM_BLOCK_SET(*block_bm_ptr,block_idx);
      // Fetch next block
      bt_block_t* const bt_block = wf_backtrace_buffer_get_block(bt_buffer,block_idx);
      const bt_block_idx_t prev_block_idx = bt_block->prev_idx;
      // Check NULL
      if (prev_block_idx != BT_BLOCK_IDX_NULL &&
          prev_block_idx >= num_compacted_blocks) {
        // Update next bt-block index and prefetch
        pf_block_idx[i] = prev_block_idx;
        BITMAP_PREFETCH_BLOCK(bitmap,prev_block_idx);
        wf_backtrace_buffer_prefetch_block(bt_buffer,prev_block_idx);
        // Next in batch
        i = (i+1) % active_blocks;
        continue;
      }
    }
    // Refill
    while (true /* !refilled */) {
      if (next_idx < num_block_idxs) {
        // Check NULL
        if (offsets[next_idx] < 0 || bt_block_idxs[next_idx] < num_compacted_blocks) {
          ++next_idx; // Next
          continue;
        }
        // Prefetch (bt-block and bt_block)
        const bt_block_idx_t block_idx = bt_block_idxs[next_idx];
        BITMAP_PREFETCH_BLOCK(bitmap,block_idx);
        wf_backtrace_buffer_prefetch_block(bt_buffer,block_idx);
        // Store
        pf_block_idx[i] = block_idx;
        ++next_idx;
        // Next in batch
        i = (i+1) % active_blocks;
        break;
      } else {
        // Take the last in the batch
        --active_blocks;
        pf_block_idx[i] = pf_block_idx[active_blocks];
        // Next in batch
        if (active_blocks) i = (i+1) % active_blocks;
        break;
      }
    }
  }
}
bt_block_idx_t wf_backtrace_buffer_compact_marked(
    wf_backtrace_buffer_t* const bt_buffer,
    bitmap_t* const bitmap,
    const int verbose) {
  // Parameters
  const int num_segments = vector_get_used(bt_buffer->segments);
  bt_block_t** const segments = vector_get_mem(bt_buffer->segments,bt_block_t*);
  const bt_block_idx_t num_compacted_blocks = bt_buffer->num_compacted_blocks;
  // Sentinels
  bt_block_idx_t read_global_pos = num_compacted_blocks;
  bt_block_idx_t write_global_pos = num_compacted_blocks;
  bt_block_idx_t read_segidx = BT_BUFFER_SEGMENT_IDX(read_global_pos);
  bt_block_idx_t read_offset = BT_BUFFER_SEGMENT_OFFSET(read_global_pos);
  bt_block_idx_t write_segidx = BT_BUFFER_SEGMENT_IDX(write_global_pos);
  bt_block_idx_t write_offset = BT_BUFFER_SEGMENT_OFFSET(write_global_pos);
  bt_block_t* read_block = segments[read_segidx] + read_offset;
  bt_block_t* write_block = segments[write_segidx] + write_offset;
  // Traverse all BT-blocks from the beginning (stored marked)
  const bt_block_idx_t max_block_idx = BT_BUFFER_IDX(bt_buffer->segment_idx,bt_buffer->segment_offset);
  while (read_global_pos < max_block_idx) {
    // Check marked block
    BITMAP_GET_BLOCK(bitmap,read_global_pos,block_bitmap_ptr);
    if (BM_BLOCK_IS_SET(*block_bitmap_ptr,read_global_pos)) {
      // Store pcigar in compacted BT-buffer
      write_block->pcigar = read_block->pcigar;
      // Translate and store index in compacted BT-buffer
      if (read_block->prev_idx == BT_BLOCK_IDX_NULL ||
          read_block->prev_idx < num_compacted_blocks) {
        write_block->prev_idx = read_block->prev_idx;
      } else {
        write_block->prev_idx = num_compacted_blocks + bitmap_erank(bitmap,read_block->prev_idx);
      }
      // Next write
      ++write_offset; ++write_block; ++write_global_pos;
      if (write_offset >= BT_BUFFER_SEGMENT_LENGTH) {
        // Next segment
        write_block = segments[++write_segidx];
        write_offset = 0;
      }
    }
    // Next read
    ++read_offset; ++read_block; ++read_global_pos;
    if (read_offset >= BT_BUFFER_SEGMENT_LENGTH) {
      // Next segment
      if (++read_segidx >= num_segments) break;
      read_block = segments[read_segidx];
      read_offset = 0;
    }
  }
  // Update next BT-buffer index
  bt_buffer->segment_offset = write_offset;
  bt_buffer->segment_idx = write_segidx;
  bt_buffer->block_next = write_block;
  bt_buffer->num_compactions++;
  // DEBUG
  if (verbose >= 3) {
    fprintf(stderr,"[WFA::BacktraceBuffer] Compacted from %lu MB to %lu MB (%2.2f%%)",
        CONVERT_B_TO_MB(read_global_pos*sizeof(bt_block_t)),
        CONVERT_B_TO_MB(write_global_pos*sizeof(bt_block_t)),
        100.0f*(float)write_global_pos/(float)read_global_pos);
  }
  // Return last index
  return write_global_pos - 1;
}



