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

#ifndef BITMAP_H_
#define BITMAP_H_

/*
 * Includes
 */
#include "../utils/commons.h"
#include "../system/mm_allocator.h"

#define BITMAP_BLOCK_ELEMENTS 64
#define BITMAP_BLOCK_MASK     0x0000000000000001ul

/*
 * Utils
 */
#define BITMAP_PREFETCH_BLOCK(bm,position) \
  PREFETCH(bm->bitmap_blocks+(position/BITMAP_BLOCK_ELEMENTS))

#define BITMAP_GET_BLOCK(bm,position,block_bitmap_ptr) \
  const uint64_t block_num = position / BITMAP_BLOCK_ELEMENTS; \
  uint64_t* const block_bitmap_ptr = &(bm->bitmap_blocks[block_num].bitmap)

#define BM_BLOCK_IS_SET(block_bitmap,position) \
  (block_bitmap & (BITMAP_BLOCK_MASK << (position % BITMAP_BLOCK_ELEMENTS)))

#define BM_BLOCK_SET(block_bitmap,position) \
  (block_bitmap |= (BITMAP_BLOCK_MASK << (position % BITMAP_BLOCK_ELEMENTS)))

/*
 * Bitmap
 */
typedef struct {
  uint64_t counter;
  uint64_t bitmap;
} bitmap_block_t;
typedef struct {
  // Bitmap
  uint64_t num_blocks;
  bitmap_block_t* bitmap_blocks;
  // MM
  mm_allocator_t* mm_allocator;
} bitmap_t;

/*
 * Setup
 */
bitmap_t* bitmap_new(
    const uint64_t length,
    mm_allocator_t* const mm_allocator);
void bitmap_delete(
    bitmap_t* const bitmap);

/*
 * Accessors
 */
void bitmap_set(
    bitmap_t* const bitmap,
    const uint64_t pos);
bool bitmap_is_set(
    bitmap_t* const bitmap,
    const uint64_t pos);
bool bitmap_check__set(
    bitmap_t* const bitmap,
    const uint64_t pos);

/*
 * Rank
 */
void bitmap_update_counters(
    bitmap_t* const bitmap);
uint64_t bitmap_erank(
    bitmap_t* const bitmap,
    const uint64_t pos);

#endif /* BITMAP_H_ */
