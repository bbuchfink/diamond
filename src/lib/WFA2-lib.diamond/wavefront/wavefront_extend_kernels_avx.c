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
 * DESCRIPTION: WaveFront-Alignment module for the "extension" of exact matches
 */

#include "wavefront_extend.h"
#include "wavefront_align.h"
#include "wavefront_compute.h"
#include "wavefront_heuristic.h"
#include "wavefront_extend_kernels.h"
#include "wavefront_extend_kernels_avx.h"

#if __AVX2__
#include <immintrin.h>
/*
 * Wavefront-Extend Inner Kernel (Scalar)
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
 * SIMD clz, use a native instruction when available (AVX512 CD or VL
 * extensions), or emulate the clz behavior.
 */
FORCE_INLINE  __m256i avx2_lzcnt_epi32(__m256i v) {
#if __AVX512CD__ && __AVX512VL__
  return _mm256_lzcnt_epi32(v);
#else
  // Emulate clz for AVX2: https://stackoverflow.com/a/58827596
  v = _mm256_andnot_si256(_mm256_srli_epi32(v,8),v); // keep 8 MSB
  v = _mm256_castps_si256(_mm256_cvtepi32_ps(v)); // convert an integer to float
  v = _mm256_srli_epi32(v,23); // shift down the exponent
  v = _mm256_subs_epu16(_mm256_set1_epi32(158),v); // undo bias
  v = _mm256_min_epi16(v,_mm256_set1_epi32(32)); // clamp at 32
  return v;
#endif
}
/*
 * Wavefront-Extend Inner Kernel (SIMD AVX2/AVX512)
 */
FORCE_NO_INLINE void wavefront_extend_matches_packed_end2end_avx2(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int lo,
    const int hi) {
  // Parameters
  wf_offset_t* const offsets = mwavefront->offsets;
  int k_min = lo;
  int k_max = hi;
  const char* pattern = wf_aligner->sequences.pattern;
  const char* text = wf_aligner->sequences.text;
  const __m256i vector_null = _mm256_set1_epi32(-1);
  const __m256i fours = _mm256_set1_epi32(4);
  const __m256i eights = _mm256_set1_epi32(8);
  const __m256i vecShuffle = _mm256_set_epi8(28,29,30,31,24,25,26,27,
                                             20,21,22,23,16,17,18,19,
                                             12,13,14,15, 8, 9,10,11,
                                             4 , 5, 6, 7, 0, 1, 2 ,3);
  const int elems_per_register = 8;
  int num_of_diagonals = k_max - k_min + 1;
  int loop_peeling_iters = num_of_diagonals % elems_per_register;
  int k;
  for (k=k_min;k<k_min+loop_peeling_iters;k++) {
    const wf_offset_t offset = offsets[k];
    if (offset < 0) continue;
    // Extend offset
    offsets[k] = wavefront_extend_matches_packed_kernel(wf_aligner,k,offset);
  }
  if (num_of_diagonals < elems_per_register) return;
  k_min += loop_peeling_iters;
  __m256i ks = _mm256_set_epi32 (
      k_min+7,k_min+6,k_min+5,k_min+4,
      k_min+3,k_min+2,k_min+1,k_min);
  // Main SIMD extension loop
  for (k=k_min;k<=k_max;k+=elems_per_register) {
    __m256i offsets_vector = _mm256_lddqu_si256 ((__m256i*)&offsets[k]);
    __m256i h_vector = offsets_vector;
    __m256i v_vector = _mm256_sub_epi32(offsets_vector,ks);
    ks =_mm256_add_epi32 (ks, eights);
    // NULL offsets will read at index 0 (avoid segfaults)
    __m256i null_mask = _mm256_cmpgt_epi32(offsets_vector,vector_null);
    v_vector = _mm256_and_si256(null_mask,v_vector);
    h_vector = _mm256_and_si256(null_mask,h_vector);
    __m256i pattern_vector = _mm256_i32gather_epi32((int const*)&pattern[0],v_vector,1);
    __m256i text_vector = _mm256_i32gather_epi32((int const*)&text[0],h_vector,1);
    // Change endianess to make the xor + clz character comparison
    pattern_vector = _mm256_shuffle_epi8(pattern_vector,vecShuffle);
    text_vector = _mm256_shuffle_epi8(text_vector,vecShuffle);
    __m256i xor_result_vector = _mm256_xor_si256(pattern_vector,text_vector);
    __m256i clz_vector = avx2_lzcnt_epi32(xor_result_vector);
    // Divide clz by 8 to get the number of equal characters
    // Assume there are sentinels on sequences so we won't count characters
    // outside the sequences
    __m256i equal_chars =  _mm256_srli_epi32(clz_vector,3);
    offsets_vector =  _mm256_add_epi32 (offsets_vector,equal_chars);
    v_vector = _mm256_add_epi32 (v_vector,fours);
    h_vector = _mm256_add_epi32 (h_vector,fours);
    // Lanes to continue == 0xffffffff, other lanes = 0
    __m256i vector_mask = _mm256_cmpeq_epi32(equal_chars,fours);
    _mm256_storeu_si256((__m256i*)&offsets[k],offsets_vector);
    int mask = _mm256_movemask_epi8(vector_mask);
    if(mask == 0) continue;
    // ctz(0) is undefined
    while (mask != 0) {
      int tz = __builtin_ctz(mask);
      int curr_k = k + (tz/4);
      const wf_offset_t offset = offsets[curr_k];
      // Extend offset
      if (offset >= 0) {
        offsets[curr_k] = wavefront_extend_matches_packed_kernel(wf_aligner,curr_k,offset);
      } else {
        offsets[curr_k] = WAVEFRONT_OFFSET_NULL;
      }
      mask &= (0xfffffff0 << tz);
    }
  }
}

#endif // AVX2
