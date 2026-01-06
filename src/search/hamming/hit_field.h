/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

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

#include <vector>
#include <cstdint>
#include <algorithm>
#if defined(_MSC_VER)
#include <intrin.h>
#endif
#include "util/math/integer.h"

using std::vector;

struct HitField {
    
    void init(size_t query_count, size_t max_target) {
        shift_ = std::max(bit_length(max_target), 8);
        word_shift_ = shift_ - 6;
        words_per_query_ = size_t(1) << word_shift_;
        data_.assign(query_count * words_per_query_, 0ULL);
    }

    void set(uint_fast32_t query, uint_fast32_t target, bool v) noexcept {
        const size_t bit = (query << shift_) | target;
        const size_t w = bit >> 6;
        const uint64_t m = 1ULL << (bit & 63);
        uint64_t& word = data_[w];
#ifdef _MSC_VER
        word ^= ((0 - (uint64_t)v) ^ word) & m;
#else
        word ^= (-(uint64_t)v ^ word) & m;
#endif       
    }

    const vector<uint_fast32_t>& hits(size_t query) {
        hits_.clear();
        const size_t base_w = query << word_shift_;
        for (size_t off = 0; off < words_per_query_; ++off) {
            uint64_t w = data_[base_w + off];
            while (w) {
                const unsigned tz = ctz64(w);
                hits_.push_back(static_cast<uint_fast32_t>((off << 6) | tz));
                w &= (w - 1);
            }
        }
        return hits_;
    }

    size_t query_count() const noexcept {
        return data_.size() >> word_shift_;
    }

private:

    static unsigned ctz64(uint64_t x) noexcept {
#if defined(_MSC_VER) && defined(_M_X64)
        unsigned long r;
        _BitScanForward64(&r, x);
        return static_cast<unsigned>(r);
#elif defined(__GNUC__) || defined(__clang__)
        return static_cast<unsigned>(__builtin_ctzll(x));
#else
        unsigned n = 0; while ((x & 1) == 0) { x >>= 1; ++n; } return n;
#endif
    }

    int shift_ = 0;
    int word_shift_ = 0;
    size_t words_per_query_ = 0;
    vector<uint64_t> data_;
    vector<uint_fast32_t> hits_;

};