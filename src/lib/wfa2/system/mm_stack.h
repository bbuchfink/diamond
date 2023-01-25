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
 * DESCRIPTION: Simple managed-memory stack that reduces memory allocation
 *   overheads. Serves memory from large memory segments and frees all memory
 *   requested at once.
 */

#ifndef MM_STACK_H_
#define MM_STACK_H_

#include "../utils/vector.h"

/*
 * Configuration
 */
#define MM_STACK_ALIGNMENT 8 // 64bits

/*
 * MM-Allocator
 */
typedef struct {
  // Memory segments
  uint64_t segment_size;          // Memory segment size (bytes)
  vector_t* segments;             // Memory segments (mm_stack_segment_t*)
  uint64_t current_segment_idx;   // Current segment being used (serving memory)
  // Malloc memory
  vector_t* malloc_requests;      // Malloc requests (void*)
  // Stack states
  vector_t* states;               // Stack saved states (mm_stack_state_t)
} mm_stack_t;

/*
 * Setup
 */
mm_stack_t* mm_stack_new(
    const uint64_t segment_size);
void mm_stack_clear(
    mm_stack_t* const mm_stack);
void mm_stack_delete(
    mm_stack_t* const mm_stack);

/*
 * Allocator
 */
void* mm_stack_allocate(
    mm_stack_t* const mm_stack,
    const uint64_t num_bytes,
    const bool zero_mem,
    const uint64_t align_bytes);

#define mm_stack_alloc(mm_stack,type) \
  ((type*)mm_stack_allocate(mm_stack,sizeof(type),false,MM_STACK_ALIGNMENT))
#define mm_stack_malloc(mm_stack,num_bytes) \
  (mm_stack_allocate(mm_stack,num_bytes,false,MM_STACK_ALIGNMENT))
#define mm_stack_calloc(mm_stack,num_elements,type,clear_mem) \
  ((type*)mm_stack_allocate(mm_stack,(num_elements)*sizeof(type),clear_mem,MM_STACK_ALIGNMENT))

#define mm_stack_uint64(mm_stack) mm_stack_malloc(mm_stack,sizeof(uint64_t))
#define mm_stack_uint32(mm_stack) mm_stack_malloc(mm_stack,sizeof(uint32_t))
#define mm_stack_uint16(mm_stack) mm_stack_malloc(mm_stack,sizeof(uint16_t))
#define mm_stack_uint8(mm_stack)  mm_stack_malloc(mm_stack,sizeof(uint8_t))

/*
 * Push/pop states
 */
void mm_stack_push(
    mm_stack_t* const mm_stack);
void mm_stack_pop(
    mm_stack_t* const mm_stack);

/*
 * Display
 */
void mm_stack_print(
    FILE* const stream,
    mm_stack_t* const mm_stack);

#endif /* MM_STACK_H_ */
