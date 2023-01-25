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

#ifndef STRING_PADDED_H_
#define STRING_PADDED_H_

/*
 * Includes
 */
#include "../utils/commons.h"
#include "../system/mm_allocator.h"

/*
 * Strings Padded
 */
typedef struct {
  // Dimensions
  int pattern_length;
  int text_length;
  // Padded strings
  char* pattern_padded;
  char* text_padded;
  // MM
  char* pattern_padded_buffer;
  char* text_padded_buffer;
  mm_allocator_t* mm_allocator;
} strings_padded_t;

/*
 * Strings (text/pattern) padded
 */
strings_padded_t* strings_padded_new(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int padding_length,
    const bool reverse_sequences,
    mm_allocator_t* const mm_allocator);
strings_padded_t* strings_padded_new_rhomb(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int padding_length,
    const bool reverse_sequences,
    mm_allocator_t* const mm_allocator);
void strings_padded_delete(
    strings_padded_t* const strings_padded);

#endif /* STRING_PADDED_H_ */
