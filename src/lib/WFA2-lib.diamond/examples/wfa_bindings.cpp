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
 * DESCRIPTION: WFA C++ Sample-Code
 */

#include <iostream>
#include <string>
#include "bindings/cpp/WFAligner.hpp"

using namespace std;
using namespace wfa;

int main(int argc,char* argv[]) {
  // Patter & Text
  string pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
  string text    = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";

  // Create a WFAligner
  WFAlignerGapAffine aligner(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);
  // Align
  aligner.alignEnd2End(pattern,text);
  cout << "WFA-Alignment returns score " << aligner.getAlignmentScore() << endl;

  // Print CIGAR
  string cigar = aligner.getAlignment();
  cout << "PATTERN: " << pattern  << endl;
  cout << "TEXT: " << text  << endl;
  cout << "CIGAR: " << cigar  << endl;
}
