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
 * DESCRIPTION: Packed CIGAR (Alignment operations in 2-bits)
 */

#ifndef WAVEFRONT_PACKED_CIGAR_H_
#define WAVEFRONT_PACKED_CIGAR_H_

#include "../utils/commons.h"
#include "wavefront_attributes.h"

/*
 * Configuration
 */
#define PCIGAR_32BITS
//#define PCIGAR_64BITS

/*
 * Packed CIGAR
 */
#define PCIGAR_NULL        0ul
#define PCIGAR_DELETION    1ul
#define PCIGAR_MISMATCH    2ul
#define PCIGAR_INSERTION   3ul

#define PCIGAR_POP_FRONT(pcigar)           pcigar = pcigar << 2
#define PCIGAR_PUSH_BACK(pcigar,operation) ((pcigar<<2) | operation)

#define PCIGAR_PUSH_BACK_INS(pcigar)   ((pcigar<<2) | PCIGAR_INSERTION)
#define PCIGAR_PUSH_BACK_DEL(pcigar)   ((pcigar<<2) | PCIGAR_DELETION)
#define PCIGAR_PUSH_BACK_MISMS(pcigar) ((pcigar<<2) | PCIGAR_MISMATCH)

#ifdef PCIGAR_32BITS
  typedef uint32_t pcigar_t;
  #define PCIGAR_MAX_LENGTH               16
  #define PCIGAR_FULL_MASK                0x40000000u /* Completely full */
  #define PCIGAR_ALMOST_FULL_MASK         0x10000000u /* 15-slots busy or more */
  #define PCIGAR_HALF_FULL_MASK           0x00010000u /*  9-slots busy or more */
  #define PCIGAR_IS_UTILISED(pcigar,mask) ((pcigar) >= mask)
  #define PCIGAR_EXTRACT(pcigar)          ((pcigar) >> 30)
  #define PCIGAR_FREE_SLOTS(pcigar)       ((pcigar)!=0) ? __builtin_clz(pcigar)/2 : PCIGAR_MAX_LENGTH;
#else
  typedef uint64_t pcigar_t;
  #define PCIGAR_MAX_LENGTH               32
  #define PCIGAR_FULL_MASK                0x4000000000000000ul /* Completely full */
  #define PCIGAR_ALMOST_FULL_MASK         0x1000000000000000ul /* 31-slots busy or more */
  #define PCIGAR_HALF_FULL_MASK           0x0000000100000000ul /* 17-slots busy or more */
  #define PCIGAR_IS_UTILISED(pcigar,mask) ((pcigar) >= mask)
  #define PCIGAR_EXTRACT(pcigar)          ((pcigar) >> 62)
  #define PCIGAR_FREE_SLOTS(pcigar)       ((pcigar)!=0) ? __builtin_clzl(pcigar)/2 : PCIGAR_MAX_LENGTH;
#endif


/*
 * Accessors
 */
int pcigar_get_length(
    const pcigar_t pcigar);
int pcigar_unpack(
    pcigar_t pcigar,
    char* cigar_buffer);

/*
 * PCIGAR unpack
 */
void pcigar_unpack_linear(
    pcigar_t pcigar,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    alignment_match_funct_t const match_funct,
    void* const match_funct_arguments,
    int* const v_pos,
    int* const h_pos,
    char* cigar_buffer,
    int* const cigar_length);
void pcigar_unpack_affine(
    pcigar_t pcigar,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    alignment_match_funct_t const match_funct,
    void* const match_funct_arguments,
    int* const v_pos,
    int* const h_pos,
    char* cigar_buffer,
    int* const cigar_length,
    affine_matrix_type* const current_matrix_type);

/*
 * Display
 */
void pcigar_print(
    FILE* const stream,
    pcigar_t pcigar);

#endif /* WAVEFRONT_PACKED_CIGAR_H_ */
