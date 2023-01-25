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

/*
 * Pragmas
 */
#ifdef __clang__
#pragma GCC diagnostic ignored "-Winitializer-overrides"
#endif

/*
 * Include
 */
#include "../utils/dna_text.h"

/*
 * Tables/Conversions Implementation
 */
const uint8_t dna_encode_table[256] =
{
  [0 ... 255] = 4,
  ['A'] = 0, ['C'] = 1, ['G'] = 2,  ['T'] = 3, ['N'] = 4,
  ['a'] = 0, ['c'] = 1, ['g'] = 2,  ['t'] = 3, ['n'] = 4,
};
const char dna_decode_table[DNA_EXTENDED_RANGE] =
{
  [ENC_DNA_CHAR_A] = DNA_CHAR_A,
  [ENC_DNA_CHAR_C] = DNA_CHAR_C,
  [ENC_DNA_CHAR_G] = DNA_CHAR_G,
  [ENC_DNA_CHAR_T] = DNA_CHAR_T,
  [ENC_DNA_CHAR_N] = DNA_CHAR_N,
};

