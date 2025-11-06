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
#include "util/math/integer.h"
#include <cstdint>
#include <algorithm>
#include <limits>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

using std::vector;

struct HitField {
    // Allocate for `query_count` queries, with per-query capacity covering targets in [0, max_target].
    // Internally rounds up to a power-of-two (#bits = 1<<shift_). Minimum block is 256 bits (shift>=8).
    void init(size_t query_count, size_t max_target) {
        shift_ = std::max(bit_length(max_target), 8);
        word_shift_ = shift_ - 6;                 // 64 == 1<<6
        words_per_query_ = size_t(1) << word_shift_;
        data_.assign(query_count * words_per_query_, 0ULL);
    }

    // Set or clear a single hit bit.
    void set(size_t query, size_t target, bool v) {
        // Caller should ensure target < (1u<<shift_). We mirror your original contract.
        const size_t bit = (query << shift_) | target;
        const size_t w = bit >> 6;
        const uint64_t m = 1ULL << (bit & 63);
        if (v)  data_[w] |= m;
        else    data_[w] &= ~m;
    }

    // Return all target indices set for the given query.
    const vector<uint_fast32_t>& hits(size_t query) {
        hits_.clear();
        const size_t base_w = query << word_shift_;

        for (size_t off = 0; off < words_per_query_; ++off) {
            uint64_t w = data_[base_w + off];
            while (w) {
                const unsigned tz = ctz64(w);             // index of lowest set bit
                hits_.push_back(static_cast<uint_fast32_t>((off << 6) | tz));
                w &= (w - 1);                              // clear lowest set bit
            }
        }
        return hits_;
    }

    size_t query_count() const {
        return data_.size() >> word_shift_;
    }

private:
    // number of bits needed to represent x (0 -> 0)
    static int bit_length(size_t x) {
        if (!x) return 0;
#if defined(_MSC_VER) && defined(_M_X64)
        unsigned long idx;
        _BitScanReverse64(&idx, static_cast<unsigned long long>(x));
        return static_cast<int>(idx) + 1;
#elif defined(__GNUC__) || defined(__clang__)
        return int(8 * sizeof(size_t) - __builtin_clzll(static_cast<unsigned long long>(x)));
#else
        // Portable fallback
        int n = 0; while (x) { ++n; x >>= 1; } return n;
#endif
    }

    static unsigned ctz64(uint64_t x) {
        // x != 0 (guaranteed by caller)
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

    int shift_ = 0;                  // bits per-query block: 1<<shift_
    int word_shift_ = 0;             // words per-query block: 1<<word_shift_ (== (1<<shift_)/64)
    size_t words_per_query_ = 0;

    vector<uint64_t> data_;          // packed bits for all queries
    vector<uint_fast32_t> hits_;     // transient buffer of hit indices for one query
};