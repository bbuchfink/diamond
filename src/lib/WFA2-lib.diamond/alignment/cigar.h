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
 * DESCRIPTION: Cigar data-structure (match/mismatch/insertion/deletion)
 */

#ifndef CIGAR_H_
#define CIGAR_H_

#include "system/mm_allocator.h"
#include "alignment/linear_penalties.h"
#include "alignment/affine_penalties.h"
#include "alignment/affine2p_penalties.h"

/*
 * CIGAR
 */
typedef struct {
  // Alignment operations
  char* operations;        // Raw alignment operations
  int max_operations;      // Maximum buffer size
  int begin_offset;        // Begin offset
  int end_offset;          // End offset
  // Score
  int score;               // Computed scored
  // CIGAR (SAM compliant)
  bool has_misms;          // Show 'X' and '=', instead of  just 'M'
  uint32_t* cigar_buffer;  // CIGAR-operations
  int cigar_length;        // Total CIGAR-operations
} cigar_t;

/*
 * Setup
 */
cigar_t* cigar_new(
    const int max_operations);
void cigar_clear(
    cigar_t* const cigar);
void cigar_resize(
    cigar_t* const cigar,
    const int max_operations);
void cigar_free(
    cigar_t* const cigar);

/*
 * Accessors
 */
bool cigar_is_null(
    cigar_t* const cigar);

int cigar_count_matches(
    cigar_t* const cigar);

void cigar_append(
    cigar_t* const cigar_dst,
    cigar_t* const cigar_src);
void cigar_append_deletion(
    cigar_t* const cigar,
    const int length);
void cigar_append_insertion(
    cigar_t* const cigar,
    const int length);

/*
 * SAM-compliant CIGAR
 */
void cigar_get_CIGAR(
    cigar_t* const cigar,
    const bool show_mismatches,
    uint32_t** const cigar_buffer,
    int* const cigar_length);

/*
 * Score
 */
int cigar_score_edit(
    cigar_t* const cigar);
int cigar_score_gap_linear(
    cigar_t* const cigar,
    linear_penalties_t* const penalties);
int cigar_score_gap_affine(
    cigar_t* const cigar,
    affine_penalties_t* const penalties);
int cigar_score_gap_affine2p(
    cigar_t* const cigar,
    affine2p_penalties_t* const penalties);

/*
 * Utils
 */
int cigar_cmp(
    cigar_t* const cigar_a,
    cigar_t* const cigar_b);
void cigar_copy(
    cigar_t* const cigar_dst,
    cigar_t* const cigar_src);

void cigar_discover_mismatches(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    cigar_t* const cigar);

void cigar_maxtrim_gap_linear(
    cigar_t* const cigar,
    linear_penalties_t* const penalties);
void cigar_maxtrim_gap_affine(
    cigar_t* const cigar,
    affine_penalties_t* const penalties);
void cigar_maxtrim_gap_affine2p(
    cigar_t* const cigar,
    affine2p_penalties_t* const penalties);

/*
 * Check
 */
bool cigar_check_alignment(
    FILE* const stream,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    cigar_t* const cigar,
    const bool verbose);

/*
 * Display
 */
void cigar_print(
    FILE* const stream,
    cigar_t* const cigar,
    const bool print_matches);
int cigar_sprint(
    char* const buffer,
    cigar_t* const cigar,
    const bool print_matches);

void cigar_print_SAM_CIGAR(
    FILE* const stream,
    cigar_t* const cigar,
    const bool show_mismatches);
int cigar_sprint_SAM_CIGAR(
    char* const buffer,
    cigar_t* const cigar,
    const bool show_mismatches);

void cigar_print_pretty(
    FILE* const stream,
    cigar_t* const cigar,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length);
int cigar_sprint_pretty(
    char* const buffer,
    cigar_t* const cigar,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length);

#endif /* CIGAR_H_ */
