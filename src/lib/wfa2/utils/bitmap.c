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
 * DESCRIPTION: Basic bitmap datastructure (static)
 */

#include "../utils/bitmap.h"
#include "../system/mm_allocator.h"

/*
 * Setup
 */
bitmap_t* bitmap_new(
    const uint64_t length,
    mm_allocator_t* const mm_allocator) {
  // Allocate
  bitmap_t* const bitmap =
      mm_allocator_alloc(mm_allocator,bitmap_t);
  bitmap->mm_allocator = mm_allocator;
  // Allocate bitmap-blocks
  const uint64_t num_blocks = DIV_CEIL(length,BITMAP_BLOCK_ELEMENTS);
  bitmap->num_blocks = num_blocks;
  bitmap->bitmap_blocks = mm_allocator_calloc(mm_allocator,num_blocks,bitmap_block_t,true);
  // Return
  return bitmap;
}
void bitmap_delete(
    bitmap_t* const bitmap) {
  // Parameters
  mm_allocator_t* const mm_allocator = bitmap->mm_allocator;
  // Free
  mm_allocator_free(mm_allocator,bitmap->bitmap_blocks);
  mm_allocator_free(mm_allocator,bitmap);
}
/*
 * Accessors
 */
void bitmap_set(
    bitmap_t* const bitmap,
    const uint64_t position) {
  // Locate block
  const uint64_t block_num = position / BITMAP_BLOCK_ELEMENTS;
  const uint64_t block_pos = position % BITMAP_BLOCK_ELEMENTS;
  // Set bitmap
  bitmap->bitmap_blocks[block_num].bitmap |= (BITMAP_BLOCK_MASK << block_pos);
}
bool bitmap_is_set(
    bitmap_t* const bitmap,
    const uint64_t position) {
  // Locate block
  const uint64_t block_num = position / BITMAP_BLOCK_ELEMENTS;
  const uint64_t block_pos = position % BITMAP_BLOCK_ELEMENTS;
  // Set bitmap
  return bitmap->bitmap_blocks[block_num].bitmap & (BITMAP_BLOCK_MASK << block_pos);
}

bool bitmap_check__set(
    bitmap_t* const bitmap,
    const uint64_t position) {
  // Locate block
  const uint64_t block_num = position / BITMAP_BLOCK_ELEMENTS;
  const uint64_t block_pos = position % BITMAP_BLOCK_ELEMENTS;
  // Check bit set
  if (bitmap->bitmap_blocks[block_num].bitmap & (BITMAP_BLOCK_MASK << block_pos)) {
    return true; // Return true (it was set)
  } else {
    // Set bitmap
    bitmap->bitmap_blocks[block_num].bitmap |= (BITMAP_BLOCK_MASK << block_pos);
    return false; // Return false (it was not set)
  }
}
/*
 * Rank
 */
void bitmap_update_counters(
    bitmap_t* const bitmap) {
  // Parameters
  const uint64_t num_blocks = bitmap->num_blocks;
  bitmap_block_t* bitmap_block = bitmap->bitmap_blocks;
  // Update all counters
  uint64_t acc_count = 0;
  uint64_t i;
  for (i=0;i<num_blocks;++i,++bitmap_block) {
    bitmap_block->counter = acc_count;
    acc_count += POPCOUNT_64(bitmap_block->bitmap);
  }
}
uint64_t bitmap_erank(
    bitmap_t* const bitmap,
    const uint64_t position) {
  // Locate block
  const uint64_t block_num = position / BITMAP_BLOCK_ELEMENTS;
  const uint64_t block_pos = position % BITMAP_BLOCK_ELEMENTS;
  // Compute e(xclusive)rank (number of bits set to one before the given position, not included)
  bitmap_block_t* const bitmap_block = bitmap->bitmap_blocks + block_num;
  const uint64_t bitmap_masked = (block_pos!=0) ? bitmap_block->bitmap << (BITMAP_BLOCK_ELEMENTS - block_pos) : 0;
  const uint64_t bitmap_count = POPCOUNT_64(bitmap_masked);
  return bitmap_block->counter + bitmap_count;
}



