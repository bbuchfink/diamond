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
 * DESCRIPTION: Padded string module to avoid handling corner conditions
 */

#include "../utils/string_padded.h"
#include "../system/mm_allocator.h"

/*
 * Strings (text/pattern) padded
 */
void strings_padded_add_padding(
    const char* const buffer,
    const int buffer_length,
    const int begin_padding_length,
    const int end_padding_length,
    const char padding_value,
    char** const buffer_padded,
    char** const buffer_padded_begin,
    const bool reverse_sequence,
    mm_allocator_t* const mm_allocator) {
  // Allocate
  const int buffer_padded_length = begin_padding_length + buffer_length + end_padding_length;
  *buffer_padded = mm_allocator_malloc(mm_allocator,buffer_padded_length);
  // Add begin padding
  memset(*buffer_padded,padding_value,begin_padding_length);
  // Copy buffer
  *buffer_padded_begin = *buffer_padded + begin_padding_length;
  if (reverse_sequence) {
    int i;
    for (i=0;i<buffer_length;i++) {
      (*buffer_padded_begin)[i] = buffer[buffer_length-1-i];
    }
  } else {
    memcpy(*buffer_padded_begin,buffer,buffer_length);
  }
  // Add end padding
  memset(*buffer_padded_begin+buffer_length,padding_value,end_padding_length);
}
strings_padded_t* strings_padded_new(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int padding_length,
    const bool reverse_sequences,
    mm_allocator_t* const mm_allocator) {
  // Allocate
  strings_padded_t* const strings_padded =
      mm_allocator_alloc(mm_allocator,strings_padded_t);
  strings_padded->mm_allocator = mm_allocator;
  // Compute padding dimensions
  const int pattern_begin_padding_length = 0;
  const int pattern_end_padding_length = padding_length;
  const int text_begin_padding_length = 0;
  const int text_end_padding_length = padding_length;
  // Add padding
  strings_padded_add_padding(
      pattern,pattern_length,
      pattern_begin_padding_length,pattern_end_padding_length,'?',
      &(strings_padded->pattern_padded_buffer),
      &(strings_padded->pattern_padded),
      reverse_sequences,mm_allocator);
  strings_padded_add_padding(
      text,text_length,
      text_begin_padding_length,text_end_padding_length,'!',
      &(strings_padded->text_padded_buffer),
      &(strings_padded->text_padded),
      reverse_sequences,mm_allocator);
  // Return
  return strings_padded;
}
strings_padded_t* strings_padded_new_rhomb(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int padding_length,
    const bool reverse_sequences,
    mm_allocator_t* const mm_allocator) {
  // Allocate
  strings_padded_t* const strings_padded =
      mm_allocator_alloc(mm_allocator,strings_padded_t);
  strings_padded->mm_allocator = mm_allocator;
  // Compute padding dimensions
  const int pattern_begin_padding_length = text_length + padding_length;
  const int pattern_end_padding_length = pattern_length + text_length + padding_length;
  const int text_begin_padding_length = padding_length;
  const int text_end_padding_length = text_length + padding_length;
  // Add padding
  strings_padded_add_padding(
      pattern,pattern_length,
      pattern_begin_padding_length,pattern_end_padding_length,'?',
      &(strings_padded->pattern_padded_buffer),
      &(strings_padded->pattern_padded),
      reverse_sequences,mm_allocator);
  strings_padded_add_padding(
      text,text_length,
      text_begin_padding_length,text_end_padding_length,'!',
      &(strings_padded->text_padded_buffer),
      &(strings_padded->text_padded),
      reverse_sequences,mm_allocator);
  // Set lengths
  strings_padded->pattern_length = pattern_length;
  strings_padded->text_length = text_length;
  // Return
  return strings_padded;
}
void strings_padded_delete(strings_padded_t* const strings_padded) {
  mm_allocator_free(strings_padded->mm_allocator,strings_padded->pattern_padded_buffer);
  mm_allocator_free(strings_padded->mm_allocator,strings_padded->text_padded_buffer);
  mm_allocator_free(strings_padded->mm_allocator,strings_padded);
}
