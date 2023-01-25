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
 * DESCRIPTION: WaveFront offset type and utils
 */

#ifndef WAVEFRONT_OFFSET_H_
#define WAVEFRONT_OFFSET_H_

#include "../utils/commons.h"

/*
 * Wavefront Offset
 */
typedef int32_t wf_offset_t;
typedef uint32_t wf_unsigned_offset_t;

/*
 * Constants
 */
#define WAVEFRONT_OFFSET_NULL (INT32_MIN/2)

/*
 * Translate k and offset to coordinates h,v
 */
#define WAVEFRONT_LENGTH(lo,hi)           ((hi)-(lo)+1)    // (lo/hi inclusive and +1 for WF[0])
#define WAVEFRONT_V(k,offset)             ((offset)-(k))
#define WAVEFRONT_H(k,offset)             (offset)
#define WAVEFRONT_ANTIDIAGONAL(k,offset)  (2*(offset)-(k))

#define DPMATRIX_DIAGONAL_NULL            INT_MAX
#define DPMATRIX_DIAGONAL(h,v)            ((h)-(v))
#define DPMATRIX_ANTIDIAGONAL(h,v)        ((h)+(v))
#define DPMATRIX_OFFSET(h,v)              (h)

#define WAVEFRONT_K_INVERSE(k,plen,tlen)  ((tlen)-(plen)-(k))

#endif /* WAVEFRONT_OFFSET_H_ */
