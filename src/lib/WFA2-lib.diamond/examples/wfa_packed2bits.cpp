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

using namespace std;
using namespace wfa;

uint8_t dna_encode_table[256];

uint8_t* stringToPacked2bits(const string& sequence) {
  const int seqLen = sequence.size();
  const int bitSeqBytes = (seqLen+7)/8;
  uint8_t* bitSeq = new uint8_t[bitSeqBytes](); // Allocate and zero
  for (int i=0;i<seqLen;++i) {
    bitSeq[i/4] |= dna_encode_table[(int)sequence[i]] << (i%4);
  }
  return bitSeq;
}
int main(int argc,char* argv[]) {
  // Patter & Text
  string pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATA";
  string text    = "TCTATACGCGCGTTTGGAGATTTAAAATAGT";

  // Create aligner
  WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);

  // Align ASCII-sequences
  aligner.alignEnd2End(pattern,text);
  cout << "WFA-Alignment (ASCII): " << aligner.getAlignmentScore() <<
      "\tCIGAR: " << aligner.getAlignment()  << endl;

  // Define encoding
  dna_encode_table['A'] = 0;
  dna_encode_table['C'] = 1;
  dna_encode_table['G'] = 2;
  dna_encode_table['T'] = 3;

  // Align packed2bits-sequences
  uint8_t* const pattern2bits = stringToPacked2bits(pattern);
  uint8_t* const text2bits = stringToPacked2bits(text);
  aligner.alignEnd2End(pattern2bits,pattern.size(),text2bits,text.size());
  cout << "WFA-Alignment (2bits): " << aligner.getAlignmentScore() <<
      "\tCIGAR: " << aligner.getAlignment()  << endl;

  // Free (because it is polite)
  delete pattern2bits;
  delete text2bits;
}
