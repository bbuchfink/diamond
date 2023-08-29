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
 * DESCRIPTION: WFA Sample-Code
 */

#include "utils/commons.h"
#include "wavefront/wavefront_align.h"

int main(int argc,char* argv[]) {
  // Patter & Text
  char* pattern = "TCTTTACTCGCGCGTTTCTTACTCGCGCGTTGGAGAAATACAATAGTGGAGAAATACAATAGTTTTTTTTTTTT";
  char* text    = "TTTTTTCTATACTGCGCGTTTTCTATACTCGCGCGTTGGAGAAATACAATAGTGGAGAAATAAAATAGT";
  // Configure alignment attributes
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine;
  attributes.affine_penalties.match = 0;
  attributes.affine_penalties.mismatch = 4;
  attributes.affine_penalties.gap_opening = 6;
  attributes.affine_penalties.gap_extension = 2;
  // Initialize Wavefront Aligner
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  // Repeat alignment (for the sake of it)
  int i;
  for (i=0;i<100000;++i) {
    // Align
    wavefront_align(wf_aligner,pattern,strlen(pattern),text,strlen(text));
    // Report
    if ((i%1000) == 0) {
      fprintf(stderr,"... done %d alignments\n",i);
    }
  }
  fprintf(stderr,"... done %d alignments\n",100000);
  // Free
  wavefront_aligner_delete(wf_aligner);
}
