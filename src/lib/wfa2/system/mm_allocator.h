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
 * VERSION: v21.02.15
 * DESCRIPTION: Simple managed-memory allocator that reduces the overhead
 *   of using malloc/calloc/free functions by allocating slabs of memory
 *   and dispatching memory segments in order.
 */

#ifndef MM_ALLOCATOR_H_
#define MM_ALLOCATOR_H_

#include "../utils/vector.h"

/*
 * Configuration
 */
//#define MM_ALLOCATOR_LOG
#define MM_ALLOCATOR_ALIGNMENT 8 // 64bits

/*
 * MM-Allocator
 */
typedef struct {
  // Metadata
  uint64_t request_ticker;        // Request ticker
  // Memory segments
  uint64_t segment_size;          // Memory segment size (bytes)
  vector_t* segments;             // Memory segments (mm_allocator_segment_t*)
  vector_t* segments_free;        // Completely free segments (mm_allocator_segment_t*)
  uint64_t current_segment_idx;   // Current segment being used (serving memory)
  // Malloc memory
  vector_t* malloc_requests;      // Malloc requests (mm_malloc_request_t)
  uint64_t malloc_requests_freed; // Total malloc request freed and still in vector
} mm_allocator_t;

/*
 * Setup
 */
mm_allocator_t* mm_allocator_new(
    const uint64_t segment_size);
void mm_allocator_clear(
    mm_allocator_t* const mm_allocator);
void mm_allocator_delete(
    mm_allocator_t* const mm_allocator);

/*
 * Allocator
 */
void* mm_allocator_allocate(
    mm_allocator_t* const mm_allocator,
    const uint64_t num_bytes,
    const bool zero_mem,
    const uint64_t align_bytes
#ifdef MM_ALLOCATOR_LOG
    ,const char* func_name,
    uint64_t line_no
#endif
    );

#ifdef MM_ALLOCATOR_LOG
#define mm_allocator_alloc(mm_allocator,type) \
  ((type*)mm_allocator_allocate(mm_allocator,sizeof(type),false,MM_ALLOCATOR_ALIGNMENT,__func__,(uint64_t)__LINE__))
#define mm_allocator_malloc(mm_allocator,num_bytes) \
  (mm_allocator_allocate(mm_allocator,num_bytes,false,MM_ALLOCATOR_ALIGNMENT,__func__,(uint64_t)__LINE__))
#define mm_allocator_calloc(mm_allocator,num_elements,type,clear_mem) \
  ((type*)mm_allocator_allocate(mm_allocator,(num_elements)*sizeof(type),clear_mem,MM_ALLOCATOR_ALIGNMENT,__func__,(uint64_t)__LINE__))
#else
#define mm_allocator_alloc(mm_allocator,type) \
  ((type*)mm_allocator_allocate(mm_allocator,sizeof(type),false,MM_ALLOCATOR_ALIGNMENT))
#define mm_allocator_malloc(mm_allocator,num_bytes) \
  (mm_allocator_allocate(mm_allocator,num_bytes,false,MM_ALLOCATOR_ALIGNMENT))
#define mm_allocator_calloc(mm_allocator,num_elements,type,clear_mem) \
  ((type*)mm_allocator_allocate(mm_allocator,(num_elements)*sizeof(type),clear_mem,MM_ALLOCATOR_ALIGNMENT))
#endif

#define mm_allocator_uint64(mm_allocator) mm_allocator_malloc(mm_allocator,sizeof(uint64_t))
#define mm_allocator_uint32(mm_allocator) mm_allocator_malloc(mm_allocator,sizeof(uint32_t))
#define mm_allocator_uint16(mm_allocator) mm_allocator_malloc(mm_allocator,sizeof(uint16_t))
#define mm_allocator_uint8(mm_allocator)  mm_allocator_malloc(mm_allocator,sizeof(uint8_t))

/*
 * Free
 */
void mm_allocator_free(
    mm_allocator_t* const mm_allocator,
    void* const memory);

/*
 * Utils
 */
void mm_allocator_get_occupation(
    mm_allocator_t* const mm_allocator,
    uint64_t* const bytes_used_malloc,
    uint64_t* const bytes_used_allocator,
    uint64_t* const bytes_free_available,
    uint64_t* const bytes_free_fragmented);

/*
 * Display
 */
void mm_allocator_print(
    FILE* const stream,
    mm_allocator_t* const mm_allocator,
    const bool display_requests);

#endif /* MM_ALLOCATOR_H_ */
