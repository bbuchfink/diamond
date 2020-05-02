#ifndef TRANSPOSE_H_
#define TRANSPOSE_H_

#include "../simd.h"

#ifdef __SSE2__
#include "transpose16x16.h"
#endif

#ifdef __AVX2__
#include "transpose32x32.h"
#endif

#endif