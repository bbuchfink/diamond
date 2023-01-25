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
 */

#ifndef SEQUENCE_BUFFER_H_
#define SEQUENCE_BUFFER_H_

#include "../utils/commons.h"
#include "../system/mm_allocator.h"

typedef struct {
  uint64_t pattern_offset;
  uint64_t pattern_length;
  uint64_t text_offset;
  uint64_t text_length;
} sequence_offset_t;
typedef struct {
  // ID
  uint64_t sequence_id;
  // Sequences
  sequence_offset_t* offsets;
  uint64_t offsets_used;
  uint64_t offsets_allocated;
  // Buffer
  char* buffer;
  uint64_t buffer_used;
  uint64_t buffer_allocated;
  // Stats
  int max_pattern_length;
  int max_text_length;
} sequence_buffer_t;

/*
 * Setup
 */
sequence_buffer_t* sequence_buffer_new(
    const uint64_t num_sequences_hint,
    const uint64_t sequence_length_hint);
void sequence_buffer_clear(
    sequence_buffer_t* const sequence_buffer);
void sequence_buffer_delete(
    sequence_buffer_t* const sequence_buffer);

/*
 * Accessors
 */
void sequence_buffer_add_pair(
    sequence_buffer_t* const sequence_buffer,
    char* const pattern,
    const uint64_t pattern_length,
    char* const text,
    const uint64_t text_length);

#endif /* SEQUENCE_BUFFER_H_ */
