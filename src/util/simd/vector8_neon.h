/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>
Arm NEON port contributed by Martin Larralde <martin.larralde@embl.de>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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
		v(vld1q_s8(p))
	{}

	operator int8x16_t() const {
		return v;
	}

	int8x16_t v;

};

template<>
struct Vector<int16_t> {

	static constexpr size_t CHANNELS = 8;

	Vector()
	{}

	Vector(const int16_t* p) :
		v(vld1q_s16(p))
	{}

	operator int16x8_t() const {
		return v;
	}

	int16x8_t v;

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