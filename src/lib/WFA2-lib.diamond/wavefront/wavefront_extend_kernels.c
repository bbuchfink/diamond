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
 * DESCRIPTION: WFA module for the "extension" of exact matches
 */

#include "wavefront_extend_kernels.h"
#include "wavefront_termination.h"

/*
 * Inner-most extend kernel (blockwise comparisons)
 */
FORCE_INLINE wf_offset_t wavefront_extend_matches_packed_kernel(
    wavefront_aligner_t* const wf_aligner,
    const int k,
    wf_offset_t offset) {
  // Fetch pattern/text blocks
  uint64_t* pattern_blocks = (uint64_t*)(wf_aligner->sequences.pattern+WAVEFRONT_V(k,offset));
  uint64_t* text_blocks = (uint64_t*)(wf_aligner->sequences.text+WAVEFRONT_H(k,offset));
  // Compare 64-bits blocks
  uint64_t cmp = *pattern_blocks ^ *text_blocks;
  while (__builtin_expect(cmp==0,0)) {
    // Increment offset (full block)
    offset += 8;
    // Next blocks
    ++pattern_blocks;
    ++text_blocks;
    // Compare
    cmp = *pattern_blocks ^ *text_blocks;
  }
  // Count equal characters
  const int equal_right_bits = __builtin_ctzl(cmp);
  const int equal_chars = DIV_FLOOR(equal_right_bits,8);
  offset += equal_chars;
  // Return extended offset
  return offset;
}
/*
 * Wavefront-Extend Inner Kernels
 *   Wavefront offset extension comparing characters
 *   Remember:
 *   - No offset is out of boundaries !(h>tlen,v>plen)
 *   - if (h==tlen,v==plen) extension won't increment (sentinels)
 */
FORCE_NO_INLINE void wavefront_extend_matches_packed_end2end(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int lo,
    const int hi) {
  wf_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=lo;k<=hi;++k) {
    // Fetch offset
    const wf_offset_t offset = offsets[k];
    if (offset == WAVEFRONT_OFFSET_NULL) continue;
    // Extend offset
    offsets[k] = wavefront_extend_matches_packed_kernel(wf_aligner,k,offset);
  }
}
FORCE_NO_INLINE wf_offset_t wavefront_extend_matches_packed_end2end_max(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int lo,
    const int hi) {
  wf_offset_t* const offsets = mwavefront->offsets;
  wf_offset_t max_antidiag = 0;
  int k;
  for (k=lo;k<=hi;++k) {
    // Fetch offset
    const wf_offset_t offset = offsets[k];
    if (offset == WAVEFRONT_OFFSET_NULL) continue;
    // Extend offset
    offsets[k] = wavefront_extend_matches_packed_kernel(wf_aligner,k,offset);
    // Compute max
    const wf_offset_t antidiag = WAVEFRONT_ANTIDIAGONAL(k,offsets[k]);
    if (max_antidiag < antidiag) max_antidiag = antidiag;
  }
  return max_antidiag;
}
FORCE_NO_INLINE bool wavefront_extend_matches_packed_endsfree(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int lo,
    const int hi) {
  // Parameters
  wf_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=lo;k<=hi;++k) {
    // Fetch offset
    wf_offset_t offset = offsets[k];
    if (offset == WAVEFRONT_OFFSET_NULL) continue;
    // Extend offset
    offset = wavefront_extend_matches_packed_kernel(wf_aligner,k,offset);
    offsets[k] = offset;
    // Check ends-free reaching boundaries
    if (wavefront_termination_endsfree(wf_aligner,mwavefront,score,k,offset)) {
      return true; // Quit (we are done)
    }
    /*
     * TODO
    const int h_pos = WAVEFRONT_H(k,offset);
    const int v_pos = WAVEFRONT_V(k,offset);
    if (h_pos >= text_length || v_pos >= pattern_length) { // FIXME Use wherever necessary
      if (wavefront_extend_endsfree_check_termination(wf_aligner,mwavefront,score,k,offset)) {
        return true; // Quit (we are done)
      }
    */
  }
  // Alignment not finished
  return false;
}
/*
 * Wavefront-Extend Inner Kernel (Custom match function)
 */
bool wavefront_extend_matches_custom(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int lo,
    const int hi,
    const bool endsfree,
    wf_offset_t* const max_antidiag) {
  // Parameters
  wavefront_sequences_t* const seqs = &wf_aligner->sequences;
  // Extend diagonally each wavefront point
  wf_offset_t* const offsets = mwavefront->offsets;
  *max_antidiag = 0;
  int k;
  for (k=lo;k<=hi;++k) {
    // Check offset
    wf_offset_t offset = offsets[k];
    if (offset == WAVEFRONT_OFFSET_NULL) continue;
    // Count equal characters
    int v = WAVEFRONT_V(k,offset);
    int h = WAVEFRONT_H(k,offset);
    while (wavefront_sequences_cmp(seqs,v,h)) {
      h++; v++; offset++;
    }
    // Update offset
    offsets[k] = offset;
    // Compute max
    const wf_offset_t antidiag = WAVEFRONT_ANTIDIAGONAL(k,offset);
    if (*max_antidiag < antidiag) *max_antidiag = antidiag;
    // Check ends-free reaching boundaries
    if (endsfree && wavefront_termination_endsfree(wf_aligner,mwavefront,score,k,offset)) {
      return true; // Quit (we are done)
    }
  }
  // Alignment not finished
  return false;
}
