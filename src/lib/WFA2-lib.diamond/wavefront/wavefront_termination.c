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
 * DESCRIPTION: WFA module to check for the termination of an alignment
 */

#include "wavefront_termination.h"

/*
 * Detect alignment termination (end of alignment)
 */
FORCE_NO_INLINE bool wavefront_termination_end2end(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int score_mod) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  const affine2p_matrix_type component_end = wf_aligner->component_end;
  const int alignment_k = DPMATRIX_DIAGONAL(text_length,pattern_length);
  const wf_offset_t alignment_offset = DPMATRIX_OFFSET(text_length,pattern_length);
  // Select end component
  switch (component_end) {
    case affine2p_matrix_M: {
      // Check diagonal/offset
      if (mwavefront->lo > alignment_k || alignment_k > mwavefront->hi) return false; // Not done
      const wf_offset_t moffset = mwavefront->offsets[alignment_k];
      if (moffset < alignment_offset) return false; // Not done
      // We are done
      wf_aligner->alignment_end_pos.score = score;
      wf_aligner->alignment_end_pos.k = alignment_k;
      wf_aligner->alignment_end_pos.offset = alignment_offset;
      return true;
    }
    case affine2p_matrix_I1: {
      // Fetch I1-wavefront & check diagonal/offset
      wavefront_t* const i1wavefront = wf_aligner->wf_components.i1wavefronts[score_mod];
      if (i1wavefront == NULL || i1wavefront->lo > alignment_k || alignment_k > i1wavefront->hi) return false; // Not done
      const wf_offset_t i1offset = i1wavefront->offsets[alignment_k];
      if (i1offset < alignment_offset) return false; // Not done
      // We are done
      wf_aligner->alignment_end_pos.score = score;
      wf_aligner->alignment_end_pos.k = alignment_k;
      wf_aligner->alignment_end_pos.offset = alignment_offset;
      return true;
    }
    case affine2p_matrix_I2: {
      // Fetch I2-wavefront & check diagonal/offset
      wavefront_t* const i2wavefront = wf_aligner->wf_components.i2wavefronts[score_mod];
      if (i2wavefront == NULL || i2wavefront->lo > alignment_k || alignment_k > i2wavefront->hi) return false; // Not done
      const wf_offset_t i2offset = i2wavefront->offsets[alignment_k];
      if (i2offset < alignment_offset) return false; // Not done
      // We are done
      wf_aligner->alignment_end_pos.score = score;
      wf_aligner->alignment_end_pos.k = alignment_k;
      wf_aligner->alignment_end_pos.offset = alignment_offset;
      return true;
    }
    case affine2p_matrix_D1: {
      // Fetch D1-wavefront & check diagonal/offset
      wavefront_t* const d1wavefront = wf_aligner->wf_components.d1wavefronts[score_mod];
      if (d1wavefront == NULL || d1wavefront->lo > alignment_k || alignment_k > d1wavefront->hi) return false; // Not done
      const wf_offset_t d1offset = d1wavefront->offsets[alignment_k];
      if (d1offset < alignment_offset) return false; // Not done
      // We are done
      wf_aligner->alignment_end_pos.score = score;
      wf_aligner->alignment_end_pos.k = alignment_k;
      wf_aligner->alignment_end_pos.offset = alignment_offset;
      return true;
    }
    case affine2p_matrix_D2: {
      // Fetch D2-wavefront & check diagonal/offset
      wavefront_t* const d2wavefront = wf_aligner->wf_components.d2wavefronts[score_mod];
      if (d2wavefront == NULL || d2wavefront->lo > alignment_k || alignment_k > d2wavefront->hi) return false; // Not done
      const wf_offset_t d2offset = d2wavefront->offsets[alignment_k];
      if (d2offset < alignment_offset) return false; // Not done
      // We are done
      wf_aligner->alignment_end_pos.score = score;
      wf_aligner->alignment_end_pos.k = alignment_k;
      wf_aligner->alignment_end_pos.offset = alignment_offset;
      return true;
    }
    default:
      break;
  }
  return false;
}
FORCE_NO_INLINE bool wavefront_termination_endsfree(
    wavefront_aligner_t* const wf_aligner,
    wavefront_t* const mwavefront,
    const int score,
    const int k,
    const wf_offset_t offset) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  // Check ends-free reaching boundaries
  const int h_pos = WAVEFRONT_H(k,offset);
  const int v_pos = WAVEFRONT_V(k,offset);
  if (h_pos >= text_length) { // Text is aligned
    // Is Pattern end-free?
    const int pattern_left = pattern_length - v_pos;
    const int pattern_end_free = wf_aligner->alignment_form.pattern_end_free;
    if (pattern_left <= pattern_end_free) {
      #ifdef WFA_PARALLEL
      #pragma omp critical
      #endif
      {
        wf_aligner->alignment_end_pos.score = score;
        wf_aligner->alignment_end_pos.k = k;
        wf_aligner->alignment_end_pos.offset = offset;
      }
      return true; // Quit (we are done)
    }
  }
  if (v_pos >= pattern_length) { // Pattern is aligned
    // Is text end-free?
    const int text_left = text_length - h_pos;
    const int text_end_free = wf_aligner->alignment_form.text_end_free;
    if (text_left <= text_end_free) {
      #ifdef WFA_PARALLEL
      #pragma omp critical
      #endif
      {
        wf_aligner->alignment_end_pos.score = score;
        wf_aligner->alignment_end_pos.k = k;
        wf_aligner->alignment_end_pos.offset = offset;
      }
      return true; // Quit (we are done)
    }
  }
  // Not done
  return false;
}
