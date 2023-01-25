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
 * DESCRIPTION: WaveFront penalties handling module
 */

#include "wavefront_penalties.h"

/*
 * Penalties adjustment
 */
void wavefront_penalties_set_indel(
    wavefront_penalties_t* const wf_penalties) {
  // Set distance model
  wf_penalties->distance_metric = indel;
  // Set penalties
  wf_penalties->match = 0;
  wf_penalties->mismatch = -1;
  wf_penalties->gap_opening1 = 1;
  wf_penalties->gap_extension1 = -1;
  wf_penalties->gap_opening2 = -1;
  wf_penalties->gap_extension2 = -1;
}
void wavefront_penalties_set_edit(
    wavefront_penalties_t* const wf_penalties) {
  // Set distance model
  wf_penalties->distance_metric = edit;
  // Set penalties
  wf_penalties->match = 0;
  wf_penalties->mismatch = 1;
  wf_penalties->gap_opening1 = 1;
  wf_penalties->gap_extension1 = -1;
  wf_penalties->gap_opening2 = -1;
  wf_penalties->gap_extension2 = -1;
}
void wavefront_penalties_set_linear(
    wavefront_penalties_t* const wf_penalties,
    linear_penalties_t* const linear_penalties) {
  // Set distance model
  wf_penalties->distance_metric = gap_linear;
  // Check base penalties
  if (linear_penalties->match > 0) {
    fprintf(stderr,"[WFA::Penalties] Match score must be negative or zero (M=%d)\n",linear_penalties->match);
    exit(1);
  } else if (linear_penalties->mismatch <= 0 || linear_penalties->indel <= 0) {
    fprintf(stderr,"[WFA::Penalties] Penalties (X=%d,D=%d,I=%d) must be (X>0,D>0,I>0)\n",
        linear_penalties->mismatch,linear_penalties->indel,linear_penalties->indel);
    exit(1);
  }
  // Set penalties (if needed, adjust using Eizenga's formula)
  if (linear_penalties->match < 0) {
    wf_penalties->match = linear_penalties->match;
    wf_penalties->mismatch = 2*linear_penalties->mismatch - 2*linear_penalties->match;
    wf_penalties->gap_opening1 = 2*linear_penalties->indel - linear_penalties->match;
  } else {
    wf_penalties->match = 0;
    wf_penalties->mismatch = linear_penalties->mismatch;
    wf_penalties->gap_opening1 = linear_penalties->indel;
  }
  // Set unused
  wf_penalties->gap_extension1 = -1;
  wf_penalties->gap_opening2 = -1;
  wf_penalties->gap_extension2 = -1;
}
void wavefront_penalties_set_affine(
    wavefront_penalties_t* const wf_penalties,
    affine_penalties_t* const affine_penalties) {
  // Set distance model
  wf_penalties->distance_metric = gap_affine;
  // Check base penalties
  if (affine_penalties->match > 0) {
    fprintf(stderr,"[WFA::Penalties] Match score must be negative or zero (M=%d)\n",affine_penalties->match);
    exit(1);
  } else if (affine_penalties->mismatch <= 0 ||
             affine_penalties->gap_opening < 0 ||
             affine_penalties->gap_extension <= 0) {
    fprintf(stderr,"[WFA::Penalties] Penalties (X=%d,O=%d,E=%d) must be (X>0,O>=0,E>0)\n",
        affine_penalties->mismatch,
        affine_penalties->gap_opening,
        affine_penalties->gap_extension);
    exit(1);
  }
  // Set penalties (if needed, adjust using Eizenga's formula)
  if (affine_penalties->match < 0) {
    wf_penalties->match = affine_penalties->match;
    wf_penalties->mismatch = 2*affine_penalties->mismatch - 2*affine_penalties->match;
    wf_penalties->gap_opening1 = 2*affine_penalties->gap_opening;
    wf_penalties->gap_extension1 = 2*affine_penalties->gap_extension - affine_penalties->match;
  } else {
    wf_penalties->match = 0;
    wf_penalties->mismatch = affine_penalties->mismatch;
    wf_penalties->gap_opening1 = affine_penalties->gap_opening;
    wf_penalties->gap_extension1 = affine_penalties->gap_extension;
  }
  // Set unused
  wf_penalties->gap_opening2 = -1;
  wf_penalties->gap_extension2 = -1;
}
void wavefront_penalties_set_affine2p(
    wavefront_penalties_t* const wf_penalties,
    affine2p_penalties_t* const affine2p_penalties) {
  // Set distance model
  wf_penalties->distance_metric = gap_affine_2p;
  // Check base penalties
  if (affine2p_penalties->match > 0) {
    fprintf(stderr,"[WFA::Penalties] Match score must be negative or zero (M=%d)\n",affine2p_penalties->match);
    exit(1);
  } else if (affine2p_penalties->mismatch <= 0 ||
             affine2p_penalties->gap_opening1 <= 0 ||
             affine2p_penalties->gap_extension1 <= 0 ||
             affine2p_penalties->gap_opening2 <= 0 ||
             affine2p_penalties->gap_extension2 <= 0) {
    fprintf(stderr,"[WFA::Penalties] Penalties (X=%d,O1=%d,E1=%d,O2=%d,E2=%d) must be (X>0,O1>=0,E1>0,O1>=0,E1>0)\n",
        affine2p_penalties->mismatch,
        affine2p_penalties->gap_opening1,
        affine2p_penalties->gap_extension1,
        affine2p_penalties->gap_opening2,
        affine2p_penalties->gap_extension2);
    exit(1);
  }
  // Set penalties (if needed, adjust using Eizenga's formula)
  if (affine2p_penalties->match < 0) {
    wf_penalties->match = affine2p_penalties->match;
    wf_penalties->mismatch = 2*affine2p_penalties->mismatch - 2*wf_penalties->match;
    wf_penalties->gap_opening1 = 2*affine2p_penalties->gap_opening1;
    wf_penalties->gap_extension1 = 2*affine2p_penalties->gap_extension1 - affine2p_penalties->match;
    wf_penalties->gap_opening2 = 2*affine2p_penalties->gap_opening2;
    wf_penalties->gap_extension2 = 2*affine2p_penalties->gap_extension2 - affine2p_penalties->match;
  } else {
    wf_penalties->match = 0;
    wf_penalties->mismatch = affine2p_penalties->mismatch;
    wf_penalties->gap_opening1 = affine2p_penalties->gap_opening1;
    wf_penalties->gap_extension1 = affine2p_penalties->gap_extension1;
    wf_penalties->gap_opening2 = affine2p_penalties->gap_opening2;
    wf_penalties->gap_extension2 = affine2p_penalties->gap_extension2;
  }
}
/*
 * Display
 */
void wavefront_penalties_print(
    FILE* const stream,
    wavefront_penalties_t* const wf_penalties) {
  // Select penalties mode
  switch (wf_penalties->distance_metric) {
    case indel:
      fprintf(stream,"(Indel)");
      break;
    case edit:
      fprintf(stream,"(Edit)");
      break;
    case gap_linear:
      fprintf(stream,"(GapLinear,%d,%d)",
          wf_penalties->mismatch,
          wf_penalties->gap_opening1);
      break;
    case gap_affine:
      fprintf(stream,"(GapAffine,%d,%d,%d)",
          wf_penalties->mismatch,
          wf_penalties->gap_opening1,
          wf_penalties->gap_extension1);
      break;
    case gap_affine_2p:
      fprintf(stream,"(GapAffine2p,%d,%d,%d,%d,%d)",
          wf_penalties->mismatch,
          wf_penalties->gap_opening1,
          wf_penalties->gap_extension1,
          wf_penalties->gap_opening2,
          wf_penalties->gap_extension2);
      break;
    default:
      break;
  }
}


