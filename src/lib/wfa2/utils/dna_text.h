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
 * DESCRIPTION: DNA text encoding/decoding utils
 */


#ifndef DNA_TEXT_H_
#define DNA_TEXT_H_

#include "../utils/commons.h"

/*
 * Range of DNA Nucleotides
 */
#define DNA_RANGE           4
#define DNA_EXTENDED_RANGE  5

#define DNA_RANGE_BITS 2

/*
 * DNA Nucleotides
 */
#define DNA_CHAR_A 'A'
#define DNA_CHAR_C 'C'
#define DNA_CHAR_G 'G'
#define DNA_CHAR_T 'T'
#define DNA_CHAR_N 'N'

/*
 * Encoded DNA Nucleotides
 */
#define ENC_DNA_CHAR_A 0
#define ENC_DNA_CHAR_C 1
#define ENC_DNA_CHAR_G 2
#define ENC_DNA_CHAR_T 3
#define ENC_DNA_CHAR_N 4

/*
 * Translation tables
 */
extern const uint8_t dna_encode_table[256];
extern const char dna_decode_table[DNA_EXTENDED_RANGE];

/*
 * Translation functions
 */
#define dna_encode(character)     (dna_encode_table[(int)(character)])
#define dna_decode(enc_char)      (dna_decode_table[(int)(enc_char)])

#endif /* DNA_TEXT_H_ */
