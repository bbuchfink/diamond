/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#pragma once
#include <stdint.h>
#include "../simd.h"

namespace DISPATCH_ARCH { namespace SIMD {

template<>
struct Vector<int8_t> {

	static constexpr size_t CHANNELS = 16;

	Vector()
	{}

	Vector(const signed char *p):
		v(_mm_loadu_si128((const __m128i*)p))
	{}

	operator __m128i() const {
		return v;
	}

	__m128i v;

};

template<>
struct Vector<int16_t> {

	static constexpr size_t CHANNELS = 8;

	Vector()
	{}

	Vector(const int16_t* p) :
		v(_mm_loadu_si128((const __m128i*)p))
	{}

	operator __m128i() const {
		return v;
	}

	__m128i v;

};

template<>
struct Vector<int32_t> {

	static constexpr size_t CHANNELS = 1;

	Vector()
	{}

	Vector(const int32_t* p) :
		v(*p)
	{}

	operator int32_t() const {
		return v;
	}

	int32_t v;

};


}}