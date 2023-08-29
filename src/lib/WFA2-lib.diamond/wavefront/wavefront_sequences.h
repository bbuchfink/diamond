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
 * DESCRIPTION: WFA module to encapsulate the input sequences
 */

#ifndef WAVEFRONT_SEQUENCES_H_
#define WAVEFRONT_SEQUENCES_H_

#include "utils/commons.h"

/*
 * Custom extend-match function, e.g.:
 *
 *   typedef struct {
 *     char* pattern;
 *     int pattern_length;
 *     char* text;
 *     int text_length;
 *   } match_function_params_t;
 *
 *   int match_function(int v,int h,void* arguments) {
 *     // Extract parameters
 *     match_function_params_t* match_arguments = (match_function_params_t*)arguments;
 *     // Check match
 *     if (v > match_arguments->pattern_length || h > match_arguments->text_length) return 0;
 *     return (match_arguments->pattern[v] == match_arguments->text[h]);
 *   }
 */
typedef int (*alignment_match_funct_t)(int,int,void*);

/*
 * Wavefront Sequences
 */
typedef enum {
  wf_sequences_ascii       = 0,
  wf_sequences_lambda      = 1,
  wf_sequences_packed2bits = 2,
} wf_sequences_mode_t;
typedef struct {
  // Mode
  wf_sequences_mode_t mode;              // Sequences mode
  bool reverse;                          // Reverse sequences
  // Current sequences & bounds
  char* pattern;                         // Pointer to current pattern sequence (padded)
  char* text;                            // Pointer to current text sequence (padded)
  int pattern_begin;                     // Pattern begin offset
  int pattern_length;                    // Pattern length
  int text_begin;                        // Text begin offset
  int text_length;                       // Text length
  // Lambda Sequence
  alignment_match_funct_t match_funct;   // Custom matching function (match(v,h,args))
  void* match_funct_arguments;           // Generic arguments passed to matching function (args)
  // Internal buffers (ASCII encoded)
  char* seq_buffer;                      // Internal buffer
  int seq_buffer_allocated;              // Internal buffer allocated
  char* pattern_buffer;                  // Source pattern sequence
  char* text_buffer;                     // Source text sequence
  int pattern_buffer_length;             // Source pattern length
  int text_buffer_length;                // Source text length
  char pattern_eos;                      // Source pattern char at EOS
  char text_eos;                         // Source pattern char at EOS
} wavefront_sequences_t;

/*
 * Setup
 */
void wavefront_sequences_allocate(
    wavefront_sequences_t* const wf_sequences);
void wavefront_sequences_free(
    wavefront_sequences_t* const wf_sequences);

/*
 * Init Sequences
 */
void wavefront_sequences_init_ascii(
    wavefront_sequences_t* const wf_sequences,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const bool reverse);
void wavefront_sequences_init_lambda(
    wavefront_sequences_t* const wf_sequences,
    alignment_match_funct_t match_funct,
    void* match_funct_arguments,
    const int pattern_length,
    const int text_length,
    const bool reverse);
void wavefront_sequences_init_packed2bits(
    wavefront_sequences_t* const wf_sequences,
    const uint8_t* const pattern,
    const int pattern_length,
    const uint8_t* const text,
    const int text_length,
    const bool reverse);


/*
 * Accessors
 */
bool wavefront_sequences_cmp(
    wavefront_sequences_t* const wf_sequences,
    const int pattern_pos,
    const int text_pos);
char wavefront_sequences_get_pattern(
    wavefront_sequences_t* const wf_sequences,
    const int position);
char wavefront_sequences_get_text(
    wavefront_sequences_t* const wf_sequences,
    const int position);

/*
 * Resize/Update
 */
void wavefront_sequences_set_bounds(
    wavefront_sequences_t* const wf_sequences,
    const int pattern_begin,
    const int pattern_end,
    const int text_begin,
    const int text_end);

#endif /* WAVEFRONT_SEQUENCES_H_ */
