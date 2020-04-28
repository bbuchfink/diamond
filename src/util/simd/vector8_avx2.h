#ifndef SIMD_VECTOR8_AVX2_H_
#define SIMD_VECTOR8_AVX2_H_

#include <stdint.h>
#include "../simd.h"

namespace SIMD {

template<>
struct Vector<int8_t> {

	Vector(const signed char* p) :
		v(_mm256_loadu_si256((const __m256i*)p))
	{}

	operator __m256i() const {
		return v;
	}

	__m256i v;

};

}

#endif