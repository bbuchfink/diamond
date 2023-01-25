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

#include "WFAligner.hpp"

extern "C" {
  #include "../../wavefront/wavefront_align.h"
}

/*
 * Namespace
 */
namespace wfa {

/*
 * General Wavefront Aligner
 */
WFAligner::WFAligner(
    const AlignmentScope alignmentScope,
    const MemoryModel memoryModel) {
  this->attributes = wavefront_aligner_attr_default;
  switch (memoryModel) {
    case MemoryHigh: this->attributes.memory_mode = wavefront_memory_high; break;
    case MemoryMed: this->attributes.memory_mode = wavefront_memory_med; break;
    case MemoryLow: this->attributes.memory_mode = wavefront_memory_low; break;
    case MemoryUltralow: this->attributes.memory_mode = wavefront_memory_ultralow; break;
    default: this->attributes.memory_mode = wavefront_memory_high; break;
  }
  this->attributes.alignment_scope = (alignmentScope==Score) ? compute_score : compute_alignment;
  //this->attributes.system.verbose = 2;
  this->wfAligner = nullptr;
}
WFAligner::~WFAligner() {
  wavefront_aligner_delete(wfAligner);
}
/*
 * Align End-to-end
 */
WFAligner::AlignmentStatus WFAligner::alignEnd2EndLambda(
    const int patternLength,
    const int textLength) {
  // Configure
  wavefront_aligner_set_alignment_end_to_end(wfAligner);
  // Align (using custom matching function)
  return (WFAligner::AlignmentStatus) wavefront_align(wfAligner,NULL,patternLength,NULL,textLength);
}
WFAligner::AlignmentStatus WFAligner::alignEnd2End(
    const char* const pattern,
    const int patternLength,
    const char* const text,
    const int textLength) {
  // Configure
  wavefront_aligner_set_alignment_end_to_end(wfAligner);
  // Align
  return (WFAligner::AlignmentStatus) wavefront_align(wfAligner,pattern,patternLength,text,textLength);
}
WFAligner::AlignmentStatus WFAligner::alignEnd2End(
    std::string& pattern,
    std::string& text) {
  // Delegate
  return alignEnd2End(pattern.c_str(),pattern.length(),text.c_str(),text.length());
}
/*
 * Align Ends-free
 */
WFAligner::AlignmentStatus WFAligner::alignEndsFreeLambda(
    const int patternLength,
    const int patternBeginFree,
    const int patternEndFree,
    const int textLength,
    const int textBeginFree,
    const int textEndFree) {
  // Configure
  wavefront_aligner_set_alignment_free_ends(wfAligner,
      patternBeginFree,patternEndFree,
      textBeginFree,textEndFree);
  // Align (using custom matching function)
  return (WFAligner::AlignmentStatus) wavefront_align(wfAligner,NULL,patternLength,NULL,textLength);
}
WFAligner::AlignmentStatus WFAligner::alignEndsFree(
    const char* const pattern,
    const int patternLength,
    const int patternBeginFree,
    const int patternEndFree,
    const char* const text,
    const int textLength,
    const int textBeginFree,
    const int textEndFree) {
  // Configure
  wavefront_aligner_set_alignment_free_ends(wfAligner,
      patternBeginFree,patternEndFree,
      textBeginFree,textEndFree);
  // Align
  return (WFAligner::AlignmentStatus) wavefront_align(wfAligner,pattern,patternLength,text,textLength);
}
WFAligner::AlignmentStatus WFAligner::alignEndsFree(
    std::string& pattern,
    const int patternBeginFree,
    const int patternEndFree,
    std::string& text,
    const int textBeginFree,
    const int textEndFree) {
  // Delegate
  return alignEndsFree(
      pattern.c_str(),pattern.length(),
      patternBeginFree,patternEndFree,
      text.c_str(),text.length(),
      textBeginFree,textEndFree);
}
/*
 * Alignment resume
 */
WFAligner::AlignmentStatus WFAligner::alignResume() {
  // Resume alignment
  return (WFAligner::AlignmentStatus) wavefront_align_resume(wfAligner);
}
/*
 * Heuristics
 */
void WFAligner::setHeuristicNone() {
  wavefront_aligner_set_heuristic_none(wfAligner);
}
void WFAligner::setHeuristicBandedStatic(
    const int band_min_k,
    const int band_max_k) {
  wavefront_aligner_set_heuristic_banded_static(
      wfAligner,band_min_k,band_max_k);
}
void WFAligner::setHeuristicBandedAdaptive(
    const int band_min_k,
    const int band_max_k,
    const int steps_between_cutoffs) {
  wavefront_aligner_set_heuristic_banded_adaptive(
      wfAligner,band_min_k,band_max_k,steps_between_cutoffs);
}
void WFAligner::setHeuristicWFadaptive(
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int steps_between_cutoffs) {
  wavefront_aligner_set_heuristic_wfadaptive(
      wfAligner,min_wavefront_length,
      max_distance_threshold,steps_between_cutoffs);
}
void WFAligner::setHeuristicXDrop(
    const int xdrop,
    const int steps_between_cutoffs) {
  wavefront_aligner_set_heuristic_xdrop(
      wfAligner,xdrop,steps_between_cutoffs);
}
void WFAligner::setHeuristicZDrop(
    const int zdrop,
    const int steps_between_cutoffs) {
  wavefront_aligner_set_heuristic_zdrop(
      wfAligner,zdrop,steps_between_cutoffs);
}
/*
 * Custom extend-match function (lambda)
 */
void WFAligner::setMatchFunct(
    int (*matchFunct)(int,int,void*),
    void* matchFunctArguments) {
  wavefront_aligner_set_match_funct(wfAligner,matchFunct,matchFunctArguments);
}
/*
 * Limits
 */
void WFAligner::setMaxAlignmentScore(
    const int maxAlignmentScore) {
  wavefront_aligner_set_max_alignment_score(
      wfAligner,maxAlignmentScore);
}
void WFAligner::setMaxMemory(
    const uint64_t maxMemoryResident,
    const uint64_t maxMemoryAbort) {
  wavefront_aligner_set_max_memory(wfAligner,maxMemoryResident,maxMemoryAbort);
}
/*
 * Accessors
 */
int WFAligner::getAlignmentScore() {
  return wfAligner->cigar->score;
}
int WFAligner::getAlignmentStatus() {
  return wfAligner->align_status.status;
}
void WFAligner::getAlignmentCigar(
    char** const cigarOperations,
    int* cigarLength) {
 *cigarOperations = wfAligner->cigar->operations + wfAligner->cigar->begin_offset;
 *cigarLength = wfAligner->cigar->end_offset - wfAligner->cigar->begin_offset;
}
std::string WFAligner::getAlignmentCigar() {
  // Fetch CIGAR
  char* buffer;
  int bufferLength;
  getAlignmentCigar(&buffer,&bufferLength);
  // Create string and return
  return std::string(buffer,bufferLength);
}
/*
 * Misc
 */
char* WFAligner::strError(
    const int wfErrorCode) {
  return wavefront_align_strerror(wfErrorCode);
}
void WFAligner::setVerbose(
    const int verbose) {
  wfAligner->system.verbose = verbose;
}
/*
 * Indel Aligner (a.k.a Longest Common Subsequence - LCS)
 */
WFAlignerIndel::WFAlignerIndel(
    const AlignmentScope alignmentScope,
    const MemoryModel memoryModel) :
        WFAligner(alignmentScope,memoryModel) {
  attributes.distance_metric = indel;
  wfAligner = wavefront_aligner_new(&attributes);
}
/*
 * Edit Aligner (a.k.a Levenshtein)
 */
WFAlignerEdit::WFAlignerEdit(
    const AlignmentScope alignmentScope,
    const MemoryModel memoryModel) :
        WFAligner(alignmentScope,memoryModel) {
  attributes.distance_metric = edit;
  wfAligner = wavefront_aligner_new(&attributes);
}
/*
 * Gap-Linear Aligner (a.k.a Needleman-Wunsch)
 */
WFAlignerGapLinear::WFAlignerGapLinear(
    const int mismatch,
    const int indel,
    const AlignmentScope alignmentScope,
    const MemoryModel memoryModel) :
        WFAligner(alignmentScope,memoryModel) {
  attributes.distance_metric = gap_linear;
  attributes.linear_penalties.match = 0;
  attributes.linear_penalties.mismatch = mismatch;
  attributes.linear_penalties.indel = indel;
  wfAligner = wavefront_aligner_new(&attributes);
}
WFAlignerGapLinear::WFAlignerGapLinear(
    const int match,
    const int mismatch,
    const int indel,
    const AlignmentScope alignmentScope,
    const MemoryModel memoryModel) :
        WFAligner(alignmentScope,memoryModel) {
  attributes.distance_metric = gap_linear;
  attributes.linear_penalties.match = match;
  attributes.linear_penalties.mismatch = mismatch;
  attributes.linear_penalties.indel = indel;
  wfAligner = wavefront_aligner_new(&attributes);
}
/*
 * Gap-Affine Aligner (a.k.a Smith-Waterman-Gotoh)
 */
WFAlignerGapAffine::WFAlignerGapAffine(
    const int mismatch,
    const int gapOpening,
    const int gapExtension,
    const AlignmentScope alignmentScope,
    const MemoryModel memoryModel) :
        WFAligner(alignmentScope,memoryModel) {
  attributes.distance_metric = gap_affine;
  attributes.affine_penalties.match = 0;
  attributes.affine_penalties.mismatch = mismatch;
  attributes.affine_penalties.gap_opening = gapOpening;
  attributes.affine_penalties.gap_extension = gapExtension;
  wfAligner = wavefront_aligner_new(&attributes);
}
WFAlignerGapAffine::WFAlignerGapAffine(
    const int match,
    const int mismatch,
    const int gapOpening,
    const int gapExtension,
    const AlignmentScope alignmentScope,
    const MemoryModel memoryModel) :
        WFAligner(alignmentScope,memoryModel) {
  attributes.distance_metric = gap_affine;
  attributes.affine_penalties.match = match;
  attributes.affine_penalties.mismatch = mismatch;
  attributes.affine_penalties.gap_opening = gapOpening;
  attributes.affine_penalties.gap_extension = gapExtension;
  wfAligner = wavefront_aligner_new(&attributes);
}
/*
 * Gap-Affine Dual-Cost Aligner (a.k.a. concave 2-pieces)
 */
WFAlignerGapAffine2Pieces::WFAlignerGapAffine2Pieces(
    const int mismatch,
    const int gapOpening1,
    const int gapExtension1,
    const int gapOpening2,
    const int gapExtension2,
    const AlignmentScope alignmentScope,
    const MemoryModel memoryModel) :
        WFAligner(alignmentScope,memoryModel) {
  attributes.distance_metric = gap_affine_2p;
  attributes.affine2p_penalties.match = 0;
  attributes.affine2p_penalties.mismatch = mismatch;
  attributes.affine2p_penalties.gap_opening1 = gapOpening1;
  attributes.affine2p_penalties.gap_extension1 = gapExtension1;
  attributes.affine2p_penalties.gap_opening2 = gapOpening2;
  attributes.affine2p_penalties.gap_extension2 = gapExtension2;
  wfAligner = wavefront_aligner_new(&attributes);
}
WFAlignerGapAffine2Pieces::WFAlignerGapAffine2Pieces(
    const int match,
    const int mismatch,
    const int gapOpening1,
    const int gapExtension1,
    const int gapOpening2,
    const int gapExtension2,
    const AlignmentScope alignmentScope,
    const MemoryModel memoryModel) :
        WFAligner(alignmentScope,memoryModel) {
  attributes.distance_metric = gap_affine_2p;
  attributes.affine2p_penalties.match = match;
  attributes.affine2p_penalties.mismatch = mismatch;
  attributes.affine2p_penalties.gap_opening1 = gapOpening1;
  attributes.affine2p_penalties.gap_extension1 = gapExtension1;
  attributes.affine2p_penalties.gap_opening2 = gapOpening2;
  attributes.affine2p_penalties.gap_extension2 = gapExtension2;
  wfAligner = wavefront_aligner_new(&attributes);
}

} /* namespace wfa */

