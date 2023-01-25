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
 * DESCRIPTION: Simple linear vector for generic type elements
 */

#include "../utils/sequence_buffer.h"

/*
 * Setup
 */
sequence_buffer_t* sequence_buffer_new(
    const uint64_t num_sequences_hint,
    const uint64_t sequence_length_hint) {
  // Alloc
  sequence_buffer_t* const sequence_buffer = malloc(sizeof(sequence_buffer_t));
  // ID
  sequence_buffer->sequence_id = 1;
  // Initialize sequences
  sequence_buffer->offsets = malloc(num_sequences_hint*sizeof(sequence_offset_t));
  sequence_buffer->offsets_used = 0;
  sequence_buffer->offsets_allocated = num_sequences_hint;
  // Initialize buffer
  const uint64_t buffer_size = num_sequences_hint*sequence_length_hint;
  sequence_buffer->buffer = malloc(buffer_size);
  sequence_buffer->buffer_used = 0;
  sequence_buffer->buffer_allocated = buffer_size;
  // Stats
  sequence_buffer->max_pattern_length = 0;
  sequence_buffer->max_text_length = 0;
  // Return
  return sequence_buffer;
}
void sequence_buffer_clear(
    sequence_buffer_t* const sequence_buffer) {
  sequence_buffer->sequence_id = 1;
  sequence_buffer->offsets_used = 0;
  sequence_buffer->buffer_used = 0;
  sequence_buffer->max_pattern_length = 0;
  sequence_buffer->max_text_length = 0;
}
void sequence_buffer_delete(
    sequence_buffer_t* const sequence_buffer) {
  free(sequence_buffer->buffer);
  free(sequence_buffer->offsets);
  free(sequence_buffer);
}
/*
 * Accessors
 */
void sequence_buffer_add_offsets(
    sequence_buffer_t* const sequence_buffer,
    const uint64_t pattern_offset,
    const uint64_t pattern_length,
    const uint64_t text_offset,
    const uint64_t text_length) {
  // Check allocated memory
  if (sequence_buffer->offsets_used == sequence_buffer->offsets_allocated) {
    const uint64_t num_offsets = (float)(sequence_buffer->offsets_used+1) * (3.0/2.0);
    sequence_buffer->offsets = realloc(sequence_buffer->offsets,num_offsets*sizeof(sequence_offset_t));
    sequence_buffer->offsets_allocated = num_offsets;
  }
  // Add sequence pair
  sequence_offset_t* const current_offsets = sequence_buffer->offsets + sequence_buffer->offsets_used;
  current_offsets->pattern_offset = pattern_offset;
  current_offsets->pattern_length = pattern_length;
  current_offsets->text_offset = text_offset;
  current_offsets->text_length = text_length;
  sequence_buffer->offsets_used += 1;
}
void sequence_buffer_add_pair(
    sequence_buffer_t* const sequence_buffer,
    char* const pattern,
    const uint64_t pattern_length,
    char* const text,
    const uint64_t text_length) {
  // Allocate memory
  const uint64_t bytes_required = pattern_length + text_length + 4; // Padding
  if (sequence_buffer->buffer_used+bytes_required > sequence_buffer->buffer_allocated) {
    const uint64_t proposed_buffer_size = (3*(sequence_buffer->buffer_used+bytes_required))/2;
    sequence_buffer->buffer_allocated = proposed_buffer_size;
    sequence_buffer->buffer = realloc(sequence_buffer->buffer,proposed_buffer_size);
  }
  // Copy pattern sequence
  char* mem_ptr = sequence_buffer->buffer + sequence_buffer->buffer_used;
  memcpy(mem_ptr,pattern,pattern_length);
  mem_ptr += pattern_length;
  mem_ptr[0] = '\0';
  mem_ptr[1] = '!';
  mem_ptr += 2; // Padding sentinels
  // Copy text sequence
  memcpy(mem_ptr,text,text_length);
  mem_ptr += text_length;
  mem_ptr[0] = '\0';
  mem_ptr[1] = '?';
  mem_ptr += 2; // Padding sentinels
  // Add sequence pair
  const uint64_t pattern_offset = sequence_buffer->buffer_used;
  const uint64_t text_offset = pattern_offset + pattern_length + 2;
  sequence_buffer_add_offsets(
      sequence_buffer,pattern_offset,
      pattern_length,text_offset,text_length);
  // Update used
  sequence_buffer->buffer_used += bytes_required;
  // Update stats
  sequence_buffer->max_pattern_length = MAX(sequence_buffer->max_pattern_length,pattern_length);
  sequence_buffer->max_text_length = MAX(sequence_buffer->max_text_length,text_length);
}

