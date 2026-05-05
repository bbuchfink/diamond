/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

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