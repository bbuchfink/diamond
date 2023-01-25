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

#include "mm_allocator.h"

/*
 * Debug
 */
//#define MM_ALLOCATOR_FORCE_MALLOC /* Delegate all requests to malloc within the mm-allocator handler */
//#define MM_ALLOCATOR_DISABLE      /* Completely disable the mm-allocator (delegate to raw malloc) */

/*
 * Constants
 */
#define MM_ALLOCATOR_SEGMENT_INITIAL_REQUESTS   10000
#define MM_ALLOCATOR_INITIAL_SEGMENTS              10
#define MM_ALLOCATOR_INITIAL_MALLOC_REQUESTS       10
#define MM_ALLOCATOR_INITIAL_STATES                10

/*
 * Allocator Segments Freed Cond
 */
#define MM_ALLOCATOR_FREED_FLAG                 0x80000000ul
#define MM_ALLOCATOR_REQUEST_IS_FREE(request)  ((request)->size & MM_ALLOCATOR_FREED_FLAG)
#define MM_ALLOCATOR_REQUEST_SET_FREE(request) ((request)->size |= MM_ALLOCATOR_FREED_FLAG)
#define MM_ALLOCATOR_REQUEST_SIZE(request)     ((request)->size & ~(MM_ALLOCATOR_FREED_FLAG))

/*
 * Reference (Header of every memory allocated)
 */
typedef struct {
  uint32_t segment_idx;
  uint32_t request_idx;
} mm_allocator_reference_t;
/*
 * Memory Request
 */
typedef struct {
  // Request
  uint32_t offset;
  uint32_t size;
  // Log
#ifdef MM_ALLOCATOR_LOG
  uint64_t timestamp;
  char* func_name;
  uint64_t line_no;
#endif
} mm_allocator_request_t;
typedef struct {
  // Request
  void* mem;
  uint64_t size;
  // Log
#ifdef MM_ALLOCATOR_LOG
  uint64_t timestamp;
  char* func_name;
  uint64_t line_no;
#endif
  // MM Reference
  mm_allocator_reference_t* reference;
} mm_malloc_request_t;
/*
 * Memory Segments
 */
typedef struct {
  // Index (ID)
  uint64_t idx;                 // Index in the segments vector
  // Memory
  uint64_t size;                // Total memory available
  void* memory;                 // Memory
  uint64_t used;                // Bytes used (offset to memory next free byte)
  // Requests
  vector_t* requests;           // Memory requests (mm_allocator_request_t)
} mm_allocator_segment_t;

/*
 * Segments
 */
mm_allocator_segment_t* mm_allocator_segment_new(
    mm_allocator_t* const mm_allocator) {
  // Allocate handler
  mm_allocator_segment_t* const segment = (mm_allocator_segment_t*) malloc(sizeof(mm_allocator_segment_t));
  // Index
  const uint64_t segment_idx = vector_get_used(mm_allocator->segments);
  segment->idx = segment_idx;
  // Memory
  segment->size = mm_allocator->segment_size;
  segment->memory = malloc(mm_allocator->segment_size);
  segment->used = 0;
  // Requests
  segment->requests = vector_new(MM_ALLOCATOR_SEGMENT_INITIAL_REQUESTS,mm_allocator_request_t);
  // Add to segments
  vector_insert(mm_allocator->segments,segment,mm_allocator_segment_t*);
  // Return
  return segment;
}
void mm_allocator_segment_clear(
    mm_allocator_segment_t* const segment) {
  segment->used = 0;
  vector_clear(segment->requests);
}
void mm_allocator_segment_delete(
    mm_allocator_segment_t* const segment) {
  vector_delete(segment->requests);
  free(segment->memory);
  free(segment);
}
mm_allocator_request_t* mm_allocator_segment_get_request(
    mm_allocator_segment_t* const segment,
    const uint64_t request_idx) {
  return vector_get_elm(segment->requests,request_idx,mm_allocator_request_t);
}
uint64_t mm_allocator_segment_get_num_requests(
    mm_allocator_segment_t* const segment) {
  return vector_get_used(segment->requests);
}
/*
 * Setup
 */
mm_allocator_t* mm_allocator_new(
    const uint64_t segment_size) {
  // Allocate handler
  mm_allocator_t* const mm_allocator = (mm_allocator_t*) malloc(sizeof(mm_allocator_t));
  mm_allocator->request_ticker = 0;
  // Segments
  mm_allocator->segment_size = segment_size;
  mm_allocator->segments = vector_new(MM_ALLOCATOR_INITIAL_SEGMENTS,mm_allocator_segment_t*);
  mm_allocator->segments_free = vector_new(MM_ALLOCATOR_INITIAL_SEGMENTS,mm_allocator_segment_t*);
  // Allocate an initial segment
#ifndef MM_ALLOCATOR_FORCE_MALLOC
#ifndef MM_ALLOCATOR_DISABLE
  mm_allocator_segment_new(mm_allocator);
#endif
#endif
  mm_allocator->current_segment_idx = 0;
  // Malloc Memory
  mm_allocator->malloc_requests = vector_new(MM_ALLOCATOR_INITIAL_MALLOC_REQUESTS,mm_malloc_request_t);
  mm_allocator->malloc_requests_freed = 0;
  // Return
  return mm_allocator;
}
void mm_allocator_clear(
    mm_allocator_t* const mm_allocator) {
  // Clear segments
  vector_clear(mm_allocator->segments_free);
  VECTOR_ITERATE(mm_allocator->segments,segment_ptr,p,mm_allocator_segment_t*) {
    mm_allocator_segment_clear(*segment_ptr); // Clear segment
    vector_insert(mm_allocator->segments_free,*segment_ptr,mm_allocator_segment_t*); // Add to free segments
  }
  mm_allocator->current_segment_idx = 0;
  // Clear malloc memory
  VECTOR_ITERATE(mm_allocator->malloc_requests,malloc_request,m,mm_malloc_request_t) {
    if (malloc_request->size > 0) free(malloc_request->mem); // Free malloc requests
  }
  vector_clear(mm_allocator->malloc_requests);
  mm_allocator->malloc_requests_freed = 0;
}
void mm_allocator_delete(
    mm_allocator_t* const mm_allocator) {
  // Free segments
  VECTOR_ITERATE(mm_allocator->segments,segment_ptr,p,mm_allocator_segment_t*) {
    mm_allocator_segment_delete(*segment_ptr);
  }
  vector_delete(mm_allocator->segments);
  vector_delete(mm_allocator->segments_free);
  // Free malloc memory
  VECTOR_ITERATE(mm_allocator->malloc_requests,malloc_request,m,mm_malloc_request_t) {
    if (malloc_request->size > 0) free(malloc_request->mem); // Free malloc requests
  }
  vector_delete(mm_allocator->malloc_requests);
  // Free handler
  free(mm_allocator);
}
/*
 * Accessors
 */
mm_allocator_segment_t* mm_allocator_get_segment(
    mm_allocator_t* const mm_allocator,
    const uint64_t segment_idx) {
  return *(vector_get_elm(mm_allocator->segments,segment_idx,mm_allocator_segment_t*));
}
mm_allocator_segment_t* mm_allocator_get_segment_free(
    mm_allocator_t* const mm_allocator,
    const uint64_t segment_idx) {
  return *(vector_get_elm(mm_allocator->segments_free,segment_idx,mm_allocator_segment_t*));
}
uint64_t mm_allocator_get_num_segments(
    mm_allocator_t* const mm_allocator) {
  return vector_get_used(mm_allocator->segments);
}
uint64_t mm_allocator_get_num_segments_free(
    mm_allocator_t* const mm_allocator) {
  return vector_get_used(mm_allocator->segments_free);
}
/*
 * Allocate
 */
mm_allocator_segment_t* mm_allocator_fetch_segment(
    mm_allocator_t* const mm_allocator,
    const uint64_t num_bytes) {
  // Fetch current segment
  mm_allocator_segment_t* const curr_segment =
      mm_allocator_get_segment(mm_allocator,mm_allocator->current_segment_idx);
  // Check overall segment size
  if (num_bytes > curr_segment->size/2) { // Never buy anything you cannot afford twice
    return NULL; // Memory request over max-request size
  }
  // Check available segment size
  if (curr_segment->used + num_bytes <= curr_segment->size) {
    return curr_segment;
  }
  // Check overall segment size
  if (num_bytes > curr_segment->size) {
    return NULL; // Memory request over segment size
  }
  // Get free segment
  const uint64_t free_segments = mm_allocator_get_num_segments_free(mm_allocator);
  if (free_segments > 0) {
    mm_allocator_segment_t* const segment =
        mm_allocator_get_segment_free(mm_allocator,free_segments-1);
    vector_dec_used(mm_allocator->segments_free);
    mm_allocator->current_segment_idx = segment->idx;
    return segment;
  }
  // Allocate new segment
  mm_allocator_segment_t* const segment = mm_allocator_segment_new(mm_allocator);
  mm_allocator->current_segment_idx = segment->idx;
  return segment;
}
void* mm_allocator_allocate(
    mm_allocator_t* const mm_allocator,
    const uint64_t num_bytes,
    const bool zero_mem,
    const uint64_t align_bytes
#ifdef MM_ALLOCATOR_LOG
    ,const char* func_name,
    uint64_t line_no
#endif
    ) {
#ifdef MM_ALLOCATOR_DISABLE
  void* memory = calloc(1,num_bytes);
  if (zero_mem) memset(memory,0,num_bytes); // Set zero
  return memory; // TODO: alignment
#else
  // Zero check
  if (num_bytes == 0) {
    fprintf(stderr,"MMAllocator error. Zero bytes requested\n");
    exit(1);
  }
  // Add payload
  const uint64_t num_bytes_allocated = num_bytes + sizeof(mm_allocator_reference_t) + align_bytes;
  // Fetch segment
#ifdef MM_ALLOCATOR_FORCE_MALLOC
  mm_allocator_segment_t* const segment = NULL; // Force malloc memory
#else
  mm_allocator_segment_t* const segment = mm_allocator_fetch_segment(mm_allocator,num_bytes_allocated);
#endif
  if (segment != NULL) {
    // Allocate memory
    void* const memory_base = segment->memory + segment->used;
    if (zero_mem) memset(memory_base,0,num_bytes_allocated); // Set zero
    // Compute aligned memory
    void* memory_aligned = memory_base + sizeof(mm_allocator_reference_t) + align_bytes;
    if (align_bytes > 0) {
      memory_aligned = memory_aligned - ((uintptr_t)memory_aligned % align_bytes);
    }
    // Set mm_reference
    mm_allocator_reference_t* const mm_reference = (mm_allocator_reference_t*)(memory_aligned - sizeof(mm_allocator_reference_t));
    mm_reference->segment_idx = segment->idx;
    mm_reference->request_idx = mm_allocator_segment_get_num_requests(segment);
    // Add request
    mm_allocator_request_t* request;
    vector_alloc_new(segment->requests,mm_allocator_request_t,request);
    request->offset = segment->used;
    request->size = num_bytes_allocated;
#ifdef MM_ALLOCATOR_LOG
    request->timestamp = (mm_allocator->request_ticker)++;
    request->func_name = (char*)func_name;
    request->line_no = line_no;
#endif
    // Update segment
    segment->used += num_bytes_allocated;
    // Return memory
    return memory_aligned;
  } else {
    // Malloc memory
    void* const memory_base = malloc(num_bytes_allocated);
    if (zero_mem) memset(memory_base,0,num_bytes_allocated); // Set zero
    // Compute aligned memory
    void* memory_aligned = memory_base + sizeof(mm_allocator_reference_t) + align_bytes;
    if (align_bytes > 0) {
      memory_aligned = memory_aligned - ((uintptr_t)memory_aligned % align_bytes);
    }
    // Set reference
    mm_allocator_reference_t* const mm_reference = (mm_allocator_reference_t*)(memory_aligned - sizeof(mm_allocator_reference_t));
    mm_reference->segment_idx = UINT32_MAX;
    mm_reference->request_idx = vector_get_used(mm_allocator->malloc_requests);
    // Add malloc-request
    mm_malloc_request_t* request;
    vector_alloc_new(mm_allocator->malloc_requests,mm_malloc_request_t,request);
    request->mem = memory_base;
    request->size = num_bytes_allocated;
#ifdef MM_ALLOCATOR_LOG
    request->timestamp = (mm_allocator->request_ticker)++;
    request->func_name = (char*)func_name;
    request->line_no = line_no;
#endif
    request->reference = mm_reference;
    // Return memory
    return memory_aligned;
  }
#endif
}
/*
 * Allocator Free
 */
void mm_allocator_free_malloc_request(
    mm_allocator_t* const mm_allocator,
    mm_allocator_reference_t* const mm_reference) {
  // Fetch request
  mm_malloc_request_t* const request =
      vector_get_elm(mm_allocator->malloc_requests,mm_reference->request_idx,mm_malloc_request_t);
  // Check double-free
  if (request->size == 0) {
    fprintf(stderr,"MMAllocator error: double free\n");
    exit(1);
  }
  // Free request
  request->size = 0;
  free(request->mem);
  ++(mm_allocator->malloc_requests_freed);
  // Check number of freed requests
  if (mm_allocator->malloc_requests_freed >= 1000) {
    // Remove freed requests
    const uint64_t num_requests = vector_get_used(mm_allocator->malloc_requests);
    mm_malloc_request_t* const requests = vector_get_mem(mm_allocator->malloc_requests,mm_malloc_request_t);
    uint64_t i, busy_requests = 0;
    for (i=0;i<num_requests;++i) {
      if (requests[i].size > 0) {
        requests[busy_requests] = requests[i];
        requests[busy_requests].reference->request_idx = busy_requests;
        ++busy_requests;
      }
    }
    vector_set_used(mm_allocator->malloc_requests,busy_requests);
    mm_allocator->malloc_requests_freed = 0;
  }
}
void mm_allocator_free_allocator_request(
    mm_allocator_t* const mm_allocator,
    mm_allocator_reference_t* const mm_reference) {
  // Fetch segment and request
  mm_allocator_segment_t* const segment =
      mm_allocator_get_segment(mm_allocator,mm_reference->segment_idx);
  mm_allocator_request_t* const request =
      mm_allocator_segment_get_request(segment,mm_reference->request_idx);
  // Check double-free
  if (MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
    fprintf(stderr,"MMAllocator error: double free\n");
    exit(1);
  }
  // Free request
  MM_ALLOCATOR_REQUEST_SET_FREE(request);
  // Free contiguous request(s) at the end of the segment
  uint64_t num_requests = mm_allocator_segment_get_num_requests(segment);
  if (mm_reference->request_idx == num_requests-1) { // Is the last request?
    --num_requests;
    mm_allocator_request_t* request =
        vector_get_mem(segment->requests,mm_allocator_request_t) + (num_requests-1);
    while (num_requests>0 && MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
      --num_requests; // Free request
      --request;
    }
    // Update segment used
    if (num_requests > 0) {
      segment->used = request->offset + request->size;
      vector_set_used(segment->requests,num_requests);
    } else {
      // Segment fully freed
      mm_allocator_segment_clear(segment); // Clear
      // Add to free segments (if it is not the current segment)
      if (segment->idx != mm_allocator->current_segment_idx) {
        vector_insert(mm_allocator->segments_free,segment,mm_allocator_segment_t*);
      }
    }
  }
}
void mm_allocator_free(
    mm_allocator_t* const mm_allocator,
    void* const memory) {
#ifdef MM_ALLOCATOR_DISABLE
  free(memory);
#else
  // Get reference
  void* const effective_memory = memory - sizeof(mm_allocator_reference_t);
  mm_allocator_reference_t* const mm_reference = (mm_allocator_reference_t*) effective_memory;
  if (mm_reference->segment_idx == UINT32_MAX) {
    // Free malloc memory
    mm_allocator_free_malloc_request(mm_allocator,mm_reference);
  } else {
    // Free allocator memory
    mm_allocator_free_allocator_request(mm_allocator,mm_reference);
  }
#endif
}
/*
 * Utils
 */
void mm_allocator_get_occupation(
    mm_allocator_t* const mm_allocator,
    uint64_t* const bytes_used_malloc,
    uint64_t* const bytes_used_allocator,
    uint64_t* const bytes_free_available,
    uint64_t* const bytes_free_fragmented) {
  // Init
  *bytes_used_malloc = 0;
  *bytes_used_allocator = 0;
  *bytes_free_available = 0;
  *bytes_free_fragmented = 0;
  // Check allocator memory
  const uint64_t num_segments = mm_allocator_get_num_segments(mm_allocator);
  int64_t segment_idx, request_idx;
  for (segment_idx=0;segment_idx<num_segments;++segment_idx) {
    mm_allocator_segment_t* const segment = mm_allocator_get_segment(mm_allocator,segment_idx);
    const uint64_t num_requests = mm_allocator_segment_get_num_requests(segment);
    bool memory_freed = true;
    for (request_idx=num_requests-1;request_idx>=0;--request_idx) {
      mm_allocator_request_t* const request = mm_allocator_segment_get_request(segment,request_idx);
      const uint64_t size = MM_ALLOCATOR_REQUEST_SIZE(request);
      if (MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
        if (memory_freed) {
          *bytes_free_available += size;
        } else {
          *bytes_free_fragmented += size;
        }
      } else {
        memory_freed = false;
        *bytes_used_allocator += size;
      }
    }
    // Account for free space at the end of the segment
    if (num_requests > 0) {
      mm_allocator_request_t* const request = mm_allocator_segment_get_request(segment,num_requests-1);
      const uint64_t bytes_free_at_end = segment->size - (request->offset+request->size);
      if (segment_idx == mm_allocator->current_segment_idx) {
        *bytes_free_available += bytes_free_at_end;
      } else {
        *bytes_free_fragmented += bytes_free_at_end;
      }
    }
  }
  // Check malloc memory
  const uint64_t num_requests = vector_get_used(mm_allocator->malloc_requests);
  mm_malloc_request_t* const requests = vector_get_mem(mm_allocator->malloc_requests,mm_malloc_request_t);
  uint64_t i;
  for (i=0;i<num_requests;++i) {
    *bytes_used_malloc += requests[i].size;
  }
}
/*
 * Display
 */
void mm_allocator_print_allocator_request(
    FILE* const stream,
    mm_allocator_request_t* const request,
    const uint64_t segment_idx,
    const uint64_t request_idx) {
      fprintf(stream,"    [#%03" PRIu64 "/%05" PRIu64 "\t%s\t@%08u\t(%" PRIu64 " Bytes)"
#ifdef MM_ALLOCATOR_LOG
          "\t%s:%" PRIu64 "\t{ts=%" PRIu64 "}"
#endif
          "\n",
          segment_idx,
          request_idx,
          MM_ALLOCATOR_REQUEST_IS_FREE(request) ? "Free]     " : "Allocated]",
          request->offset,
          (uint64_t)MM_ALLOCATOR_REQUEST_SIZE(request)
#ifdef MM_ALLOCATOR_LOG
          ,request->func_name,
          request->line_no,
          request->timestamp
#endif
      );
}
void mm_allocator_print_malloc_request(
    FILE* const stream,
    mm_malloc_request_t* const request) {
      fprintf(stream,"    [@%p" PRIu64 "\t(%" PRIu64 " Bytes)"
#ifdef MM_ALLOCATOR_LOG
          "\t%s:%" PRIu64 "\t{ts=%" PRIu64 "}"
#endif
          "\n",
          request->mem,
          request->size
#ifdef MM_ALLOCATOR_LOG
          ,request->func_name,
          request->line_no,
          request->timestamp
#endif
      );
}
void mm_allocator_print_allocator_requests(
    FILE* const stream,
    mm_allocator_t* const mm_allocator,
    const bool compact_free) {
  // Print allocator memory
  uint64_t segment_idx, request_idx;
  uint64_t free_block = 0;
  bool has_requests = false;
  fprintf(stream,"  => MMAllocator.requests\n");
  const uint64_t num_segments = mm_allocator_get_num_segments(mm_allocator);
  for (segment_idx=0;segment_idx<num_segments;++segment_idx) {
    mm_allocator_segment_t* const segment = mm_allocator_get_segment(mm_allocator,segment_idx);
    const uint64_t num_requests = mm_allocator_segment_get_num_requests(segment);
    for (request_idx=0;request_idx<num_requests;++request_idx) {
      mm_allocator_request_t* const request = mm_allocator_segment_get_request(segment,request_idx);
      if (compact_free) {
        if (MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
          free_block += MM_ALLOCATOR_REQUEST_SIZE(request);
        } else {
          if (free_block > 0) {
            fprintf(stream,"    [n/a\tFree]      \t(%" PRIu64 " Bytes)\n",free_block);
            free_block = 0;
          }
          mm_allocator_print_allocator_request(stream,request,segment_idx,request_idx);
          has_requests = true;
        }
      } else {
        mm_allocator_print_allocator_request(stream,request,segment_idx,request_idx);
        has_requests = true;
      }
    }
  }
  if (!has_requests) {
    fprintf(stream,"    -- No requests --\n");
  }
  // Print malloc memory
  fprintf(stream,"  => MMMalloc.requests\n");
  const uint64_t num_requests = vector_get_used(mm_allocator->malloc_requests);
  mm_malloc_request_t* const requests = vector_get_mem(mm_allocator->malloc_requests,mm_malloc_request_t);
  uint64_t i;
  for (i=0;i<num_requests;++i) {
    if (requests[i].size > 0) {
      mm_allocator_print_malloc_request(stream,requests+i);
    }
  }
  if (num_requests == 0) {
    fprintf(stream,"    -- No requests --\n");
  }
}
void mm_allocator_print(
    FILE* const stream,
    mm_allocator_t* const mm_allocator,
    const bool display_requests) {
  // Print header
  fprintf(stream,"MMAllocator.report\n");
  // Print segment information
  const uint64_t num_segments = mm_allocator_get_num_segments(mm_allocator);
  const uint64_t segment_size = mm_allocator->segment_size;
  fprintf(stream,"  => Segments.allocated %" PRIu64 "\n",num_segments);
  fprintf(stream,"  => Segments.size      %" PRIu64 " MB\n",segment_size/(1024*1024));
  fprintf(stream,"  => Memory.available   %" PRIu64 " MB\n",num_segments*(segment_size/(1024*1024)));
  // Print memory information
  uint64_t bytes_used_malloc, bytes_used_allocator;
  uint64_t bytes_free_available, bytes_free_fragmented;
  mm_allocator_get_occupation(mm_allocator,&bytes_used_malloc,&bytes_used_allocator,&bytes_free_available,&bytes_free_fragmented);
  const float bytes_total = num_segments * segment_size;
  const uint64_t bytes_free = bytes_free_available + bytes_free_fragmented;
  fprintf(stream,"    => Memory.used   %" PRIu64 " (%2.1f %%)\n",
      bytes_used_allocator,100.0f*(float)bytes_used_allocator/bytes_total);
  fprintf(stream,"    => Memory.free   %" PRIu64 " (%2.1f %%)\n",
      bytes_free,100.0f*(float)bytes_free/bytes_total);
  fprintf(stream,"      => Memory.free.available  %" PRIu64 " (%2.1f %%)\n",
      bytes_free_available,100.0f*(float)bytes_free_available/bytes_total);
  fprintf(stream,"      => Memory.free.fragmented %" PRIu64 " (%2.1f %%)\n",
      bytes_free_fragmented,100.0f*(float)bytes_free_fragmented/bytes_total);
  fprintf(stream,"    => Memory.malloc %" PRIu64 "\n",bytes_used_malloc);
  // Print memory requests
  if (display_requests) {
    mm_allocator_print_allocator_requests(stream,mm_allocator,false);
  }
}



