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

#include "mm_stack.h"

/*
 * Debug
 */
//#define MM_STACK_FORCE_MALLOC

/*
 * Constants
 */
#define MM_STACK_INITIAL_SEGMENTS              10
#define MM_STACK_INITIAL_MALLOC_REQUESTS       10
#define MM_STACK_INITIAL_STATES                10

/*
 * Stack state
 */
typedef struct {
  uint64_t segment_idx;
  uint64_t segment_used;
  uint64_t num_malloc_requests;
} mm_stack_state_t;

/*
 * Memory Segments
 */
typedef struct {
  uint64_t size;                // Total memory available
  void* memory;                 // Memory
  uint64_t used;                // Bytes used (offset to memory next free byte)
} mm_stack_segment_t;

/*
 * Segments
 */
mm_stack_segment_t* mm_stack_segment_new(
    mm_stack_t* const mm_stack) {
  // Allocate handler
  mm_stack_segment_t* const segment = (mm_stack_segment_t*) malloc(sizeof(mm_stack_segment_t));
  // Memory
  segment->size = mm_stack->segment_size;
  segment->memory = malloc(mm_stack->segment_size);
  segment->used = 0;
  // Add to segments
  vector_insert(mm_stack->segments,segment,mm_stack_segment_t*);
  // Return
  return segment;
}
void mm_stack_segment_clear(
    mm_stack_segment_t* const segment) {
  segment->used = 0;
}
void mm_stack_segment_delete(
    mm_stack_segment_t* const segment) {
  free(segment->memory);
  free(segment);
}
/*
 * Setup
 */
mm_stack_t* mm_stack_new(
    const uint64_t segment_size) {
  // Allocate handler
  mm_stack_t* const mm_stack = (mm_stack_t*) malloc(sizeof(mm_stack_t));
  // Memory segments
  mm_stack->segments = vector_new(MM_STACK_INITIAL_SEGMENTS,mm_stack_segment_t*);
  mm_stack->segment_size = segment_size;
#ifndef MM_STACK_FORCE_MALLOC
  mm_stack_segment_new(mm_stack);
#endif
  mm_stack->current_segment_idx = 0;
  // Malloc memory
  mm_stack->malloc_requests = vector_new(MM_STACK_INITIAL_MALLOC_REQUESTS,void*);
  // Stack states
  mm_stack->states = vector_new(MM_STACK_INITIAL_STATES,mm_stack_state_t);
  // Return
  return mm_stack;
}
void mm_stack_clear(
    mm_stack_t* const mm_stack) {
  // Clear first memory segment and discard the rest
  mm_stack_segment_t* const segment = *vector_get_elm(mm_stack->segments,0,mm_stack_segment_t*);
  mm_stack_segment_clear(segment);
  mm_stack->current_segment_idx = 0;
  // Free malloc memory
  VECTOR_ITERATE(mm_stack->malloc_requests,mem_ptr,m,void*) {
    free(*mem_ptr);
  }
  vector_clear(mm_stack->malloc_requests);
  // Clear states
  vector_clear(mm_stack->states);
}
void mm_stack_delete(
    mm_stack_t* const mm_stack) {
  // Delete memory segments
  VECTOR_ITERATE(mm_stack->segments,segment_ptr,p,mm_stack_segment_t*) {
    mm_stack_segment_delete(*segment_ptr);
  }
  vector_delete(mm_stack->segments);
  // Free malloc memory
  VECTOR_ITERATE(mm_stack->malloc_requests,mem_ptr,m,void*) {
    free(*mem_ptr);
  }
  vector_delete(mm_stack->malloc_requests);
  // Clear states
  vector_delete(mm_stack->states);
  // Free handler
  free(mm_stack);
}
/*
 * Allocator
 */
mm_stack_segment_t* mm_stack_fetch_segment(
    mm_stack_t* const mm_stack,
    const uint64_t num_bytes) {
  // Fetch current segment
  mm_stack_segment_t* const curr_segment =
      *vector_get_elm(mm_stack->segments,mm_stack->current_segment_idx,mm_stack_segment_t*);
//  // Check overall segment size
//  if (num_bytes > curr_segment->size/2) { // Never buy anything you cannot afford twice
//    return NULL; // Memory request over max-request size
//  }
  // Check available segment size
  if (curr_segment->used + num_bytes <= curr_segment->size) {
    return curr_segment;
  }
  // Check overall segment size
  if (num_bytes > curr_segment->size) {
    return NULL; // Memory request over segment size
  }
  // Get free segment
  const uint64_t num_segments = vector_get_used(mm_stack->segments);
  ++(mm_stack->current_segment_idx);
  if (mm_stack->current_segment_idx < num_segments) {
    // Get next segment
    mm_stack_segment_t* const segment =
        *vector_get_elm(mm_stack->segments,mm_stack->current_segment_idx,mm_stack_segment_t*);
    // Clear
    mm_stack_segment_clear(segment);
    // Return
    return segment;
  }
  // Add new segment
  return mm_stack_segment_new(mm_stack);
}
void* mm_stack_allocate(
    mm_stack_t* const mm_stack,
    const uint64_t num_bytes,
    const bool zero_mem,
    const uint64_t align_bytes) {
  // Zero check
  if (num_bytes == 0) {
    fprintf(stderr,"MMStack error. Zero bytes requested\n");
    exit(1);
  }
  // Add payload
  const uint64_t num_bytes_allocated = num_bytes + align_bytes;
  // Fetch segment
#ifdef MM_STACK_FORCE_MALLOC
  mm_stack_segment_t* const segment = NULL; // Force malloc memory
#else
  mm_stack_segment_t* const segment = mm_stack_fetch_segment(mm_stack,num_bytes_allocated);
#endif
  // Allocate memory
  void* memory_base ;
  if (segment != NULL) {
    // Segment-memory
    memory_base = segment->memory + segment->used;
    if (zero_mem) memset(memory_base,0,num_bytes_allocated); // Set zero
    segment->used += num_bytes_allocated; // Update segment
  } else {
    // Malloc-memory
    memory_base = malloc(num_bytes_allocated);
    if (zero_mem) memset(memory_base,0,num_bytes_allocated); // Set zero
    // Add malloc-request
    vector_insert(mm_stack->malloc_requests,memory_base,void*);
  }
  // Check alignment
  if (align_bytes == 0) return memory_base;
  // Align memory request
  void* memory_aligned = memory_base + align_bytes;
  memory_aligned = memory_aligned - ((uintptr_t)memory_aligned % align_bytes);
  return memory_aligned;
}
/*
 * Push/pop states
 */
void mm_stack_push(
    mm_stack_t* const mm_stack) {
  // Get new stack-state
  mm_stack_state_t* stack_state;
  vector_alloc_new(mm_stack->states,mm_stack_state_t,stack_state);
  // Store current state
  mm_stack_segment_t* const current_segment =
      *vector_get_elm(mm_stack->segments,mm_stack->current_segment_idx,mm_stack_segment_t*);
  stack_state->segment_idx = mm_stack->current_segment_idx;
  stack_state->segment_used = current_segment->used;
  stack_state->num_malloc_requests = vector_get_used(mm_stack->malloc_requests);
}
void mm_stack_pop(
    mm_stack_t* const mm_stack) {
  // Get last stack-state
  mm_stack_state_t* const stack_state = vector_get_last_elm(mm_stack->states,mm_stack_state_t);
  vector_dec_used(mm_stack->states);
  // Restore segment-memory state
  mm_stack->current_segment_idx = stack_state->segment_idx;
  mm_stack_segment_t* const current_segment =
      *(vector_get_elm(mm_stack->segments,stack_state->segment_idx,mm_stack_segment_t*));
  current_segment->used = stack_state->segment_used;
  // Restore malloc-memory state (free requests)
  const uint64_t total_malloc_requests = vector_get_used(mm_stack->malloc_requests);
  void** const malloc_requests = vector_get_mem(mm_stack->malloc_requests,void*);
  uint64_t i;
  for (i=stack_state->num_malloc_requests;i<total_malloc_requests;++i) {
    free(*(malloc_requests+i)); // Free
  }
  vector_set_used(mm_stack->malloc_requests,stack_state->num_malloc_requests);
}
/*
 * Display
 */
void mm_stack_print(
    FILE* const stream,
    mm_stack_t* const mm_stack) {
  // Print header
  fprintf(stream,"MMStack.report\n");
  // Print segment information
  const uint64_t num_segments = vector_get_used(mm_stack->segments);
  const uint64_t segment_size = mm_stack->segment_size;
  fprintf(stream,"  => Segments.allocated %" PRIu64 "\n",num_segments);
  fprintf(stream,"  => Segments.size      %" PRIu64 " MB\n",segment_size/(1024*1024));
  fprintf(stream,"  => Memory.available   %" PRIu64 " MB\n",num_segments*(segment_size/(1024*1024)));
}


