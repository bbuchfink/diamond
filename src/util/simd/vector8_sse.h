#ifndef SIMD_VECTOR8_SSE_H_
#define SIMD_VECTOR8_SSE_H_

#include <stdint.h>
#include "../simd.h"

namespace SIMD {

template<>
struct Vector<int8_t> {

	Vector(const signed char *p):
		v(_mm_loadu_si128((const __m128i*)p))
	{}

	operator __m128i() const {
		return v;
	}

	__m128i v;

};

}

#endif