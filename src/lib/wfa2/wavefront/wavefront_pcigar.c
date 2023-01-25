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
 * DESCRIPTION: Packed CIGAR (Alignment operations in 2-bits)
 */

#include "wavefront_pcigar.h"

/*
 * Packed CIGAR operations
 */
typedef struct {
  char operation;
  int inc_v;
  int inc_h;
  affine_matrix_type matrix_type;
} pcigar_op_t;
pcigar_op_t pcigar_lut[4] = {
    { .operation = '?', .inc_v = 0, .inc_h = 0, .matrix_type = affine_matrix_M }, // 00 - None
    { .operation = 'D', .inc_v = 1, .inc_h = 0, .matrix_type = affine_matrix_D }, // 01 - DELETION
    { .operation = 'X', .inc_v = 1, .inc_h = 1, .matrix_type = affine_matrix_M }, // 10 - MISMATCH
    { .operation = 'I', .inc_v = 0, .inc_h = 1, .matrix_type = affine_matrix_I }, // 11 - INSERTION
};
// Precomputed string of Matches
char matches_lut[8] = "MMMMMMMM";
#define CIGAR_8MATCHES_UINT64 *((uint64_t*)matches_lut)

/*
 * Accessors
 */
int pcigar_get_length(
    const pcigar_t pcigar) {
  int cigar_length = PCIGAR_MAX_LENGTH;
  if (!PCIGAR_IS_UTILISED(pcigar,PCIGAR_FULL_MASK)) {
    const int free_slots = PCIGAR_FREE_SLOTS(pcigar);
    cigar_length -= free_slots;
  }
  return cigar_length;
}
int pcigar_unpack(
    pcigar_t pcigar,
    char* cigar_buffer) {
  // Compute pcigar length and shift to the end of the word
  int pcigar_length = PCIGAR_MAX_LENGTH;
  if (!PCIGAR_IS_UTILISED(pcigar,PCIGAR_FULL_MASK)) {
    const int free_slots = PCIGAR_FREE_SLOTS(pcigar);
    pcigar_length -= free_slots;
    pcigar <<= free_slots*2;
  }
  // Unpack BT-blocks
  int i;
  for (i=0;i<pcigar_length;++i) {
    // Extract next CIGAR operation
    const int cigar_op = (int)PCIGAR_EXTRACT(pcigar); // Extract
    PCIGAR_POP_FRONT(pcigar); // Shift
    // Add operation using LUT
    pcigar_op_t* const op = pcigar_lut + cigar_op;
    *(cigar_buffer++) = op->operation;
  }
  return pcigar_length;
}
/*
 * PCIGAR extend exact-matches
 */
int pcigar_unpack_extend(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    int v,
    int h,
    char* cigar_buffer) {
  int num_matches = 0;
  // Fetch pattern/text blocks
  uint64_t* pattern_blocks = (uint64_t*)(pattern+v);
  uint64_t* text_blocks = (uint64_t*)(text+h);
  uint64_t pattern_block = *pattern_blocks;
  uint64_t text_block = *text_blocks;
  // Compare 64-bits blocks
  uint64_t cmp = pattern_block ^ text_block;
  while (cmp==0 && (v+8) < pattern_length && (h+8) < text_length) {
    // Increment offset
    v += 8;
    h += 8;
    num_matches += 8;
    // Dump matches in block
    *((uint64_t*)cigar_buffer) = CIGAR_8MATCHES_UINT64;
    cigar_buffer += 8;
    // Next blocks
    ++pattern_blocks;
    ++text_blocks;
    // Fetch & Compare
    pattern_block = *pattern_blocks;
    text_block = *text_blocks;
    cmp = pattern_block ^ text_block;
  }
  // Count equal characters
  num_matches += __builtin_ctzl(cmp)/8;
  *((uint64_t*)cigar_buffer) = CIGAR_8MATCHES_UINT64;
  // Return total matches
  return num_matches;
}
int pcigar_unpack_extend_custom(
    const int pattern_length,
    const int text_length,
    alignment_match_funct_t const match_funct,
    void* const match_funct_arguments,
    int v,
    int h,
    char* cigar_buffer) {
  int num_matches = 0;
  while (v < pattern_length && h < text_length) {
    // Check match
    if (!match_funct(v,h,match_funct_arguments)) break;
    ++v; ++h;
    // Increment matches
    *cigar_buffer = 'M';
    ++cigar_buffer;
    ++num_matches;
  }
  return num_matches;
}
/*
 * PCIGAR unpack
 */
void pcigar_unpack_linear(
    pcigar_t pcigar,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    alignment_match_funct_t const match_funct,
    void* const match_funct_arguments,
    int* const v_pos,
    int* const h_pos,
    char* cigar_buffer,
    int* const cigar_length) {
  // Parameters
  char* const cigar_buffer_base = cigar_buffer;
  // Compute pcigar length and shift to the end of the word
  int pcigar_length = PCIGAR_MAX_LENGTH;
  if (!PCIGAR_IS_UTILISED(pcigar,PCIGAR_FULL_MASK)) {
    const int free_slots = PCIGAR_FREE_SLOTS(pcigar);
    pcigar_length -= free_slots;
    pcigar <<= free_slots*2;
  }
  // Unpack BT-blocks
  int v = *v_pos, h = *h_pos, i;
  for (i=0;i<pcigar_length;++i) {
    // Extend exact-matches
    int num_matches;
    if (match_funct != NULL) { // Custom extend-match function
      num_matches = pcigar_unpack_extend_custom(
          pattern_length,text_length,
          match_funct,match_funct_arguments,v,h,cigar_buffer);
    } else {
      num_matches = pcigar_unpack_extend(
          pattern,pattern_length,text,text_length,v,h,cigar_buffer);
    }
    // Update location
    v += num_matches;
    h += num_matches;
    cigar_buffer += num_matches;
    // Extract next CIGAR operation
    const int cigar_op = (int)PCIGAR_EXTRACT(pcigar); // Extract
    PCIGAR_POP_FRONT(pcigar); // Shift
    // Add operation using LUT
    pcigar_op_t* const op = pcigar_lut + cigar_op;
    *(cigar_buffer++) = op->operation;
    v += op->inc_v;
    h += op->inc_h;
  }
  // Update length/positions
  *cigar_length = cigar_buffer - cigar_buffer_base;
  *v_pos = v;
  *h_pos = h;
}
void pcigar_unpack_affine(
    pcigar_t pcigar,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    alignment_match_funct_t const match_funct,
    void* const match_funct_arguments,
    int* const v_pos,
    int* const h_pos,
    char* cigar_buffer,
    int* const cigar_length,
    affine_matrix_type* const current_matrix_type) {
  // Parameters
  char* const cigar_buffer_base = cigar_buffer;
  // Compute pcigar length and shift to the end of the word
  int pcigar_length = PCIGAR_MAX_LENGTH;
  if (!PCIGAR_IS_UTILISED(pcigar,PCIGAR_FULL_MASK)) {
    const int free_slots = PCIGAR_FREE_SLOTS(pcigar);
    pcigar_length -= free_slots;
    pcigar <<= free_slots*2;
  }
  // Unpack BT-blocks
  affine_matrix_type matrix_type = *current_matrix_type;
  int v = *v_pos, h = *h_pos, i;
  for (i=0;i<pcigar_length;++i) {
    // Extend exact-matches
    if (matrix_type == affine_matrix_M) { // Extend only on the M-wavefront
      int num_matches;
      if (match_funct != NULL) { // Custom extend-match function
        num_matches = pcigar_unpack_extend_custom(
            pattern_length,text_length,
            match_funct,match_funct_arguments,v,h,cigar_buffer);
      } else {
        num_matches = pcigar_unpack_extend(
            pattern,pattern_length,text,text_length,v,h,cigar_buffer);
      }
      // Update location
      v += num_matches;
      h += num_matches;
      cigar_buffer += num_matches;
    }
    // Extract next CIGAR operation
    const int cigar_op = (int)PCIGAR_EXTRACT(pcigar); // Extract
    PCIGAR_POP_FRONT(pcigar); // Shift
    // Add operation using LUT
    pcigar_op_t* const op = pcigar_lut + cigar_op;
    // Avoid X after I/D (used to encode gap-close)
    if (matrix_type != affine_matrix_M && op->operation=='X') {
      matrix_type = affine_matrix_M;
      continue;
    }
    // Add operation
    *(cigar_buffer++) = op->operation;
    v += op->inc_v;
    h += op->inc_h;
    matrix_type = op->matrix_type;
  }
  // Update length/positions
  *cigar_length = cigar_buffer - cigar_buffer_base;
  *v_pos = v;
  *h_pos = h;
  *current_matrix_type = matrix_type;
}
/*
 * Display
 */
void pcigar_print(
    FILE* const stream,
    pcigar_t pcigar) {
  char cigar_buffer[64];
  const int pcigar_length = pcigar_unpack(pcigar,cigar_buffer);
  cigar_buffer[pcigar_length] = '\0';
  fprintf(stream,"%s",cigar_buffer);
}
