#ifndef SIMD_VECTOR_H_
#define SIMD_VECTOR_H_

#include "../simd.h"

#ifdef __AVX2__
#include "vector8_avx2.h"
#elif defined(__SSE2__)
#include "vector8_sse.h"
#endif

#endif