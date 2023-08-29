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
 * DESCRIPTION: WFA C++ Sample-Code. Using 2bit-encoded sequences
 */

#include <iostream>
#include <string>

#include "bindings/cpp/WFAligner.hpp"
#include "external/ksw2.h"

using namespace std;
using namespace wfa;

void align_ksw2(
    const char *tseq,
    const char *qseq,
    int sc_mch,
    int sc_mis,
    int gapo,
    int gape) {
    int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
  int8_t mat[25] = {
      static_cast<int8_t>(a),static_cast<int8_t>(b),static_cast<int8_t>(b),static_cast<int8_t>(b),0,
      static_cast<int8_t>(b),static_cast<int8_t>(a),static_cast<int8_t>(b),static_cast<int8_t>(b),0,
      static_cast<int8_t>(b),static_cast<int8_t>(b),static_cast<int8_t>(a),static_cast<int8_t>(b),0,
      static_cast<int8_t>(b),static_cast<int8_t>(b),static_cast<int8_t>(b),static_cast<int8_t>(a),0,
      0,0,0,0,0};
  int tl = strlen(tseq), ql = strlen(qseq);
  uint8_t *ts, *qs, c[256];
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  memset(c, 4, 256);
  c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
  c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
  ts = (uint8_t*)malloc(tl);
  qs = (uint8_t*)malloc(ql);
  for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
  for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];

  ksw_extz2_sse(nullptr, ql, qs, tl, ts, 5, mat, gapo, gape, -1, 30, 100, 0x40 , &ez);
  for (i = 0; i < ez.n_cigar; ++i) {
      printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
  }
  putchar('\n');
  free(ez.cigar); free(ts); free(qs);
}

void align_wfa(
    const char *tseq,
    const char *qseq,
    int sc_mch,
    int sc_mis,
    int gapo,
    int gape) {
  // Parameters
  int tl = strlen(tseq), ql = strlen(qseq);

  // Create aligner
  WFAlignerGapAffine aligner(-sc_mch,sc_mis,gapo,gape,WFAligner::Alignment);
  aligner.setHeuristicNone();
  aligner.setHeuristicZDrop(30,1); //20

  // Align sequences
  aligner.alignExtension(tseq,tl,qseq,ql);
  cout << aligner.getCIGARString(false) << "\t(" << aligner.getAlignmentScore() << ")" << endl;
}

int main(int argc,char* argv[]) {
  // Penalties
  int match = 3, mismatch = 3, gapopen = 4, gapextend = 1;

  // This is just a very short example database for the later examples
  const char* target =
      "TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGT"
      "CGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTTGTCCGGGTG"
      "TGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTAC"
      "GTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCC"
      "TCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATT"
      "CAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATA"
      "AAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAA"
      "CACTAAACATAGCAGTGGTGTTACCCGTGAACTCATGCGTGAGCTTAACGGAGGGGCATACACTCGCTAT";

  /*
   * Example 1: The query is a perfect subpart of the target, the alignment should end right after the query ends
   */
  const char* query_perfect_match =
      "TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGT"
      "CGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGC";
  std::cout << "[Example-1] Perfect Match: " << endl;
  std::cout << "  KSW\t"; align_ksw2(target,query_perfect_match,match,mismatch,gapopen,gapextend);
  std::cout << "  WFA\t"; align_wfa(target,query_perfect_match,match,mismatch,gapopen,gapextend);

  /*
   * Example 2: The query contains 37 Insertions, the alignment should continue over this gap,
   *   as there are not enough gaps to trigger the z-drop
   */
  const char* query_insertions =
      "TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGT"
      "CGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCAAAAAAAAAAACGCGCGCGCGCGCCAAAAAAAAGCGCAGCTTACGGTTTCGTCCGTGTTG"
      "CAGCCGATCATCAGCACATCTAGGTTTTGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAAC";
  std::cout << "[Example-2] Query Insertions: " << endl;
  std::cout << "  KSW\t"; align_ksw2(target,query_insertions,match,mismatch,gapopen,gapextend);
  std::cout << "  WFA\t"; align_wfa(target,query_insertions,match,mismatch,gapopen,gapextend);

  /*
   * Example 3: The query contains a long run of Insertions but only 10 Matches before the insertions,
   *   this time the heuristic should be dropped and the alignment should stop after the matches
   */
  const char* query_insertions_long = "TTGTAGATCTAGGGGGGGGCACAGCCTACGCATACATCCCCCCCCCCAAAAAAAAGGGGGGGGGGAAA"
      "AAATTTTTTGGGGGGGGAAAAAACCCGCGCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAAC";
  std::cout << "[Example-3] Query long  Insertion: " << endl;
  std::cout << "  KSW\t"; align_ksw2(target,query_insertions_long,match,mismatch,gapopen,gapextend);
  std::cout << "  WFA\t"; align_wfa(target,query_insertions_long,match,mismatch,gapopen,gapextend);
}
