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
 * DESCRIPTION: C++ bindings for the WaveFront Alignment modules
 */

#ifndef BINDINGS_CPP_WFALIGNER_HPP_
#define BINDINGS_CPP_WFALIGNER_HPP_

#include <string>

extern "C" {
  #include "../../wavefront/wavefront_aligner.h"
}

/*
 * Namespace
 */
namespace wfa {

/*
 * General Wavefront Aligner
 */
class WFAligner {
public:
  // Configuration
  enum MemoryModel {
    MemoryHigh,
    MemoryMed,
    MemoryLow,
    MemoryUltralow,
  };
  enum AlignmentScope {
    Score,
    Alignment,
  };
  enum AlignmentStatus {
    StatusSuccessful = WF_STATUS_SUCCESSFUL,
    StatusUnfeasible = WF_STATUS_UNFEASIBLE,
    StatusMaxScoreReached = WF_STATUS_MAX_SCORE_REACHED,
    StatusOOM = WF_STATUS_OOM,
  };
  // Align End-to-end
  AlignmentStatus alignEnd2EndLambda(
      const int patternLength,
      const int textLength);
  AlignmentStatus alignEnd2End(
      const char* const pattern,
      const int patternLength,
      const char* const text,
      const int textLength);
  AlignmentStatus alignEnd2End(
      std::string& pattern,
      std::string& text);
  // Align Ends-free
  AlignmentStatus alignEndsFreeLambda(
      const int patternLength,
      const int patternBeginFree,
      const int patternEndFree,
      const int textLength,
      const int textBeginFree,
      const int textEndFree);
  AlignmentStatus alignEndsFree(
      const char* const pattern,
      const int patternLength,
      const int patternBeginFree,
      const int patternEndFree,
      const char* const text,
      const int textLength,
      const int textBeginFree,
      const int textEndFree);
  AlignmentStatus alignEndsFree(
      std::string& pattern,
      const int patternBeginFree,
      const int patternEndFree,
      std::string& text,
      const int textBeginFree,
      const int textEndFree);
  // Alignment resume
  AlignmentStatus alignResume();
  // Heuristics
  void setHeuristicNone();
  void setHeuristicBandedStatic(
      const int band_min_k,
      const int band_max_k);
  void setHeuristicBandedAdaptive(
      const int band_min_k,
      const int band_max_k,
      const int steps_between_cutoffs = 1);
  void setHeuristicWFadaptive(
      const int min_wavefront_length,
      const int max_distance_threshold,
      const int steps_between_cutoffs = 1);
  void setHeuristicXDrop(
      const int xdrop,
      const int steps_between_cutoffs = 1);
  void setHeuristicZDrop(
      const int zdrop,
      const int steps_between_cutoffs = 1);
  // Custom extend-match function (lambda)
  void setMatchFunct(
      int (*matchFunct)(int,int,void*),
      void* matchFunctArguments);
  // Limits
  void setMaxAlignmentScore(
      const int maxAlignmentScore);
  void setMaxMemory(
      const uint64_t maxMemoryResident,
      const uint64_t maxMemoryAbort);
  // Accessors
  int getAlignmentScore();
  int getAlignmentStatus();
  void getAlignmentCigar(
      char** const cigarOperations,
      int* cigarLength);
  std::string getAlignmentCigar();
  // Misc
  char* strError(
      const int wfErrorCode);
  void setVerbose(
      const int verbose);
protected:
  wavefront_aligner_attr_t attributes;
  wavefront_aligner_t* wfAligner;
  // Setup
  WFAligner(
      const AlignmentScope alignmentScope,
      const MemoryModel memoryModel = MemoryHigh);
  ~WFAligner();
private:
  WFAligner(const WFAligner&);
};
/*
 * Indel Aligner (a.k.a Longest Common Subsequence - LCS)
 */
class WFAlignerIndel : public WFAligner {
public:
  WFAlignerIndel(
      const AlignmentScope alignmentScope,
      const MemoryModel memoryModel = MemoryHigh);
};
/*
 * Edit Aligner (a.k.a Levenshtein)
 */
class WFAlignerEdit : public WFAligner {
public:
  WFAlignerEdit(
      const AlignmentScope alignmentScope,
      const MemoryModel memoryModel = MemoryHigh);
};
/*
 * Gap-Linear Aligner (a.k.a Needleman-Wunsch)
 */
class WFAlignerGapLinear : public WFAligner {
public:
  WFAlignerGapLinear(
      const int mismatch,
      const int indel,
      const AlignmentScope alignmentScope,
      const MemoryModel memoryModel = MemoryHigh);
  WFAlignerGapLinear(
      const int match,
      const int mismatch,
      const int indel,
      const AlignmentScope alignmentScope,
      const MemoryModel memoryModel = MemoryHigh);
};
/*
 * Gap-Affine Aligner (a.k.a Smith-Waterman-Gotoh)
 */
class WFAlignerGapAffine : public WFAligner {
public:
  WFAlignerGapAffine(
      const int mismatch,
      const int gapOpening,
      const int gapExtension,
      const AlignmentScope alignmentScope,
      const MemoryModel memoryModel = MemoryHigh);
  WFAlignerGapAffine(
      const int match,
      const int mismatch,
      const int gapOpening,
      const int gapExtension,
      const AlignmentScope alignmentScope,
      const MemoryModel memoryModel = MemoryHigh);
};
/*
 * Gap-Affine Dual-Cost Aligner (a.k.a. concave 2-pieces)
 */
class WFAlignerGapAffine2Pieces : public WFAligner {
public:
  WFAlignerGapAffine2Pieces(
      const int mismatch,
      const int gapOpening1,
      const int gapExtension1,
      const int gapOpening2,
      const int gapExtension2,
      const AlignmentScope alignmentScope,
      const MemoryModel memoryModel = MemoryHigh);
  WFAlignerGapAffine2Pieces(
      const int match,
      const int mismatch,
      const int gapOpening1,
      const int gapExtension1,
      const int gapOpening2,
      const int gapExtension2,
      const AlignmentScope alignmentScope,
      const MemoryModel memoryModel = MemoryHigh);
};

} /* namespace wfa */

#endif /* BINDINGS_CPP_WFALIGNER_HPP_ */
