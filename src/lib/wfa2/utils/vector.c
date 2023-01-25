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
 * VERSION: v20.08.25
 * DESCRIPTION: Simple linear vector (generic type elements)
 */

#include "vector.h"

/*
 * Constants
 */
#define VECTOR_EXPAND_FACTOR (3.0/2.0)

/*
 * Setup
 */
vector_t* vector_new_(
    const uint64_t num_initial_elements,
    const uint64_t element_size) {
  vector_t* const vector_buffer = (vector_t*) malloc(sizeof(vector_t));
  vector_buffer->element_size = element_size;
  vector_buffer->elements_allocated = num_initial_elements;
  vector_buffer->memory = malloc(num_initial_elements*element_size);
  if (!vector_buffer->memory) {
    fprintf(stderr,"Could not create new vector (%" PRIu64 " bytes requested)",
        num_initial_elements*element_size);
    exit(1);
  }
  vector_buffer->used = 0;
  return vector_buffer;
}
void vector_delete(
    vector_t* const vector) {
  free(vector->memory);
  free(vector);
}
void vector_cast(
    vector_t* const vector,
    const uint64_t element_size) {
  vector->elements_allocated = (vector->elements_allocated*vector->element_size)/element_size;
  vector->element_size = element_size;
  vector->used = 0;
}
void vector_reserve(
    vector_t* const vector,
    const uint64_t num_elements,
    const bool zero_mem) {
  if (vector->elements_allocated < num_elements) {
    const uint64_t proposed = (float)vector->elements_allocated*VECTOR_EXPAND_FACTOR;
    vector->elements_allocated = num_elements>proposed?num_elements:proposed;
    vector->memory = realloc(vector->memory,vector->elements_allocated*vector->element_size);
    if (!vector->memory) {
      fprintf(stderr,"Could not reserve vector (%" PRIu64 " bytes requested)",
          vector->elements_allocated*vector->element_size);
      exit(1);
    }
  }
  if (zero_mem) {
    memset(vector->memory+vector->used*vector->element_size,0,
        (vector->elements_allocated-vector->used)*vector->element_size);
  }
}
/*
 * Accessors
 */
#ifdef VECTOR_DEBUG
void* vector_get_mem_element(
    vector_t* const vector,
    const uint64_t position,
    const uint64_t element_size) {
  if (position >= (vector)->used) {
    fprintf(stderr,"Vector position out-of-range [0,%"PRIu64")",(vector)->used);
    exit(1);
  }
  return vector->memory + (position*element_size);
}
#endif
/*
 * Miscellaneous
 */
void vector_copy(
    vector_t* const vector_to,
    vector_t* const vector_from) {
  // Prepare
  vector_cast(vector_to,vector_from->element_size);
  vector_reserve(vector_to,vector_from->used,false);
  // Copy
  vector_set_used(vector_to,vector_from->used);
  memcpy(vector_to->memory,vector_from->memory,vector_from->used*vector_from->element_size);
}
vector_t* vector_dup(
    vector_t* const vector_src) {
  vector_t* const vector_cpy = vector_new_(vector_src->used,vector_src->element_size);
  // Copy
  vector_set_used(vector_cpy,vector_src->used);
  memcpy(vector_cpy->memory,vector_src->memory,vector_src->used*vector_src->element_size);
  return vector_cpy;
}

