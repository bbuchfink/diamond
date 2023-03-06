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

#include "wavefront_sequences.h"

/*
 * Configuration
 */
#define WF_SEQUENCES_PADDING     64
#define WF_SEQUENCES_PATTERN_EOS '!'
#define WF_SEQUENCES_TEXT_EOS    '?'

/*
 * Setup
 */
void wavefront_sequences_allocate(
    wavefront_sequences_t* const wf_sequences) {
  // Mode
  wf_sequences->mode = wf_sequences_ascii;
  wf_sequences->reverse = false;
  // Source sequences
  wf_sequences->seq_buffer = NULL;
  wf_sequences->seq_buffer_allocated = 0;
  // Current state
  wf_sequences->pattern = NULL;
  wf_sequences->text = NULL;
}
void wavefront_sequences_free(
    wavefront_sequences_t* const wf_sequences) {
  // Free internal buffers
  if (wf_sequences->seq_buffer != NULL) free(wf_sequences->seq_buffer);
}
/*
 * Init Sequences
 */
void wavefront_sequences_init_allocate(
    wavefront_sequences_t* const wf_sequences,
    const int pattern_length,
    const int text_length) {
  // Compute dimensions
  const int buffer_size = pattern_length + text_length + 3*WF_SEQUENCES_PADDING;
  // Check internal buffer allocated
  if (wf_sequences->seq_buffer_allocated < buffer_size) {
    // Free
    if (wf_sequences->seq_buffer != NULL) free(wf_sequences->seq_buffer);
    // Allocate
    const int proposed_size = buffer_size + buffer_size/2;
    wf_sequences->seq_buffer = malloc(proposed_size);
    wf_sequences->seq_buffer_allocated = proposed_size;
  }
  // Assign memory
  wf_sequences->pattern_buffer = wf_sequences->seq_buffer + WF_SEQUENCES_PADDING;
  wf_sequences->text_buffer = wf_sequences->seq_buffer + WF_SEQUENCES_PADDING + pattern_length + WF_SEQUENCES_PADDING;
}
void wavefront_sequences_init_copy(
    char* const buffer_dst,
    const char* const sequence,
    const int sequence_length,
    const int padding_length,
    const char padding_value,
    const bool reverse) {
  // Copy sequence
  if (reverse) {
    int i;
    for (i=0;i<sequence_length;i++) {
      buffer_dst[i] = sequence[sequence_length-1-i];
    }
  } else {
    memcpy(buffer_dst,sequence,sequence_length);
  }
  // Add end padding
  buffer_dst[sequence_length] = padding_value;
}
void wavefront_sequences_init_decode2bits(
    char* const buffer_dst,
    const uint8_t* const sequence,
    const int sequence_length,
    const int padding_length,
    const char padding_value,
    const bool reverse) {
  // Parameters
  const char dna_packed2bits_decode[4] = {'A','C','G','T'};
  // Compute dimensions
  const int num_words = DIV_CEIL(sequence_length,8);
  int buffer_pos = (reverse) ? sequence_length-1 : 0;
  // Decode and copy sequence
  int word_num;
  for (word_num=0;word_num<num_words;++word_num) {
    // Fetch next word
    const uint8_t word = sequence[word_num];
    // Decode 4 letters packed
    const char letter0 = dna_packed2bits_decode[(word    & 3)];
    const char letter1 = dna_packed2bits_decode[(word>>2 & 3)];
    const char letter2 = dna_packed2bits_decode[(word>>4 & 3)];
    const char letter3 = dna_packed2bits_decode[(word>>6 & 3)];
    if (reverse) {
      buffer_dst[buffer_pos  ] = letter0;
      buffer_dst[buffer_pos-1] = letter1;
      buffer_dst[buffer_pos-2] = letter2;
      buffer_dst[buffer_pos-3] = letter3;
      buffer_pos -= 4;
    } else {
      buffer_dst[buffer_pos  ] = letter0;
      buffer_dst[buffer_pos+1] = letter1;
      buffer_dst[buffer_pos+2] = letter2;
      buffer_dst[buffer_pos+3] = letter3;
      buffer_pos += 4;
    }
  }
  // Add end padding
  buffer_dst[sequence_length] = padding_value;
}
void wavefront_sequences_init_ascii(
    wavefront_sequences_t* const wf_sequences,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const bool reverse) {
  // Mode
  wf_sequences->mode = wf_sequences_ascii;
  wf_sequences->reverse = reverse;
  // Allocate buffers
  wavefront_sequences_init_allocate(wf_sequences,pattern_length,text_length);
  // Copy internal sequences
  wavefront_sequences_init_copy(wf_sequences->pattern_buffer,
      pattern,pattern_length,WF_SEQUENCES_PADDING,WF_SEQUENCES_PATTERN_EOS,reverse);
  wf_sequences->pattern_buffer_length = pattern_length;
  wavefront_sequences_init_copy(wf_sequences->text_buffer,
      text,text_length,WF_SEQUENCES_PADDING,WF_SEQUENCES_TEXT_EOS,reverse);
  wf_sequences->text_buffer_length = text_length;
  // Set pattern
  wf_sequences->pattern = wf_sequences->pattern_buffer;
  wf_sequences->pattern_begin = 0;
  wf_sequences->pattern_length = pattern_length;
  wf_sequences->pattern_eos = wf_sequences->pattern[pattern_length];
  // Set text
  wf_sequences->text = wf_sequences->text_buffer;
  wf_sequences->text_begin = 0;
  wf_sequences->text_length = text_length;
  wf_sequences->text_eos = wf_sequences->text[text_length];
}
void wavefront_sequences_init_lambda(
    wavefront_sequences_t* const wf_sequences,
    alignment_match_funct_t match_funct,
    void* match_funct_arguments,
    const int pattern_length,
    const int text_length,
    const bool reverse) {
  // Mode
  wf_sequences->mode = wf_sequences_lambda;
  wf_sequences->reverse = reverse;
  // Set sequences' length
  wf_sequences->pattern = NULL;
  wf_sequences->text = NULL;
  wf_sequences->pattern_begin = 0;
  wf_sequences->pattern_length = pattern_length;
  wf_sequences->text_begin = 0;
  wf_sequences->text_length = text_length;
  // Internals
  wf_sequences->match_funct = match_funct;
  wf_sequences->match_funct_arguments = match_funct_arguments;
}
void wavefront_sequences_init_packed2bits(
    wavefront_sequences_t* const wf_sequences,
    const uint8_t* const pattern,
    const int pattern_length,
    const uint8_t* const text,
    const int text_length,
    const bool reverse) {
  // Mode
  wf_sequences->mode = wf_sequences_ascii;
  wf_sequences->reverse = reverse;
  // Allocate buffers
  wavefront_sequences_init_allocate(wf_sequences,pattern_length,text_length);
  // Copy internal sequences
  wavefront_sequences_init_decode2bits(wf_sequences->pattern_buffer,
      pattern,pattern_length,WF_SEQUENCES_PADDING,WF_SEQUENCES_PATTERN_EOS,reverse);
  wf_sequences->pattern_buffer_length = pattern_length;
  wavefront_sequences_init_decode2bits(wf_sequences->text_buffer,
      text,text_length,WF_SEQUENCES_PADDING,WF_SEQUENCES_TEXT_EOS,reverse);
  wf_sequences->text_buffer_length = text_length;
  // Set pattern
  wf_sequences->pattern = wf_sequences->pattern_buffer;
  wf_sequences->pattern_begin = 0;
  wf_sequences->pattern_length = pattern_length;
  wf_sequences->pattern_eos = wf_sequences->pattern[pattern_length];
  // Set text
  wf_sequences->text = wf_sequences->text_buffer;
  wf_sequences->text_begin = 0;
  wf_sequences->text_length = text_length;
  wf_sequences->text_eos = wf_sequences->text[text_length];
}
/*
 * Accessors
 */
bool wavefront_sequences_cmp(
    wavefront_sequences_t* const wf_sequences,
    const int pattern_pos,
    const int text_pos) {
  // Select mode
  if (wf_sequences->mode == wf_sequences_lambda) {
    // Custom function to compare sequences
    alignment_match_funct_t match_funct = wf_sequences->match_funct;
    void* match_funct_arguments = wf_sequences->match_funct_arguments;
    // Check coordinates (EOS)
    const int pattern_length = wf_sequences->pattern_length;
    const int text_length = wf_sequences->text_length;
    if (pattern_pos >= pattern_length || text_pos >= text_length) return false;
    // Compare using lambda (given coordinates)
    const int pattern_begin = wf_sequences->pattern_begin;
    const int text_begin = wf_sequences->text_begin;
    if (wf_sequences->reverse) {
      const int pattern_end = pattern_begin + pattern_length - 1;
      const int text_end = text_begin + text_length - 1;
      return match_funct(pattern_end-pattern_pos,text_end-text_pos,match_funct_arguments);
    } else {
      return match_funct(pattern_begin+pattern_pos,text_begin+text_pos,match_funct_arguments);
    }
  } else {
    // Compare regular strings
    return wf_sequences->pattern[pattern_pos] == wf_sequences->text[text_pos];
  }
}
char wavefront_sequences_get_pattern(
    wavefront_sequences_t* const wf_sequences,
    const int position) {
  if (wf_sequences->mode == wf_sequences_lambda) {
    return '-';
  } else {
    return wf_sequences->pattern[position];
  }
}
char wavefront_sequences_get_text(
    wavefront_sequences_t* const wf_sequences,
    const int position) {
  if (wf_sequences->mode == wf_sequences_lambda) {
    return '-';
  } else {
    return wf_sequences->text[position];
  }
}
/*
 * Resize/Update
 */
void wavefront_sequences_set_bounds(
    wavefront_sequences_t* const wf_sequences,
    const int pattern_begin,
    const int pattern_end,
    const int text_begin,
    const int text_end) {
  // Select mode
  if (wf_sequences->mode != wf_sequences_lambda) {
    // Restore previous EOS char
    const int pattern_length_old = wf_sequences->pattern_length;
    const int text_length_old = wf_sequences->text_length;
    wf_sequences->pattern[pattern_length_old] = wf_sequences->pattern_eos;
    wf_sequences->text[text_length_old] = wf_sequences->text_eos;
    // Focus on the new section of the sequences
    if (wf_sequences->reverse) {
      // Compare given coordinates
      wf_sequences->pattern = wf_sequences->pattern_buffer + (wf_sequences->pattern_buffer_length - pattern_end);
      wf_sequences->text = wf_sequences->text_buffer + (wf_sequences->text_buffer_length - text_end);
    } else {
      wf_sequences->pattern = wf_sequences->pattern_buffer + pattern_begin;
      wf_sequences->text = wf_sequences->text_buffer + text_begin;
    }
    // Save EOS char and truncate sequence
    const int pattern_length_new = pattern_end - pattern_begin;
    const int text_length_new = text_end - text_begin;
    wf_sequences->pattern_eos = wf_sequences->pattern[pattern_length_new];
    wf_sequences->text_eos = wf_sequences->text[text_length_new];
    wf_sequences->pattern[pattern_length_new] = WF_SEQUENCES_PATTERN_EOS;
    wf_sequences->text[text_length_new] = WF_SEQUENCES_TEXT_EOS;
  }
  // Set bounds
  wf_sequences->pattern_begin = pattern_begin;
  wf_sequences->pattern_length = pattern_end - pattern_begin;
  wf_sequences->text_begin = text_begin;
  wf_sequences->text_length = text_end - text_begin;
}

