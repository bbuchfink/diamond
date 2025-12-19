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
#include <cmath>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include "../hash_function.h"
#include "../intrin.h"

struct HyperLogLog {

    HyperLogLog(int precision = 10) : p(precision), m(1 << precision), registers(m, 0) {
        if (p < 4 || p > 20) throw std::invalid_argument("Precision must be between 4 and 20");
        compute_alpha();
    }

    void add(int64_t x) {
        uint64_t hash = MurmurHash()(x);
        process_hash(hash);
    }

    int64_t estimate() const {
        double sum = 0.0;
        int zeros = 0;
        bool all_zero = true;

        for (uint8_t r : registers) {
            sum += 1.0 / (1ULL << r);
            if (r != 0) all_zero = false;
            if (r == 0) zeros++;
        }

        if (all_zero) return 0.0;

        double z = 1.0 / sum;
        double e = alpha * m * m * z;

        if (e <= 2.5 * m && zeros > 0) {
            e = m * std::log(static_cast<double>(m) / zeros);
        }

        return (int64_t)std::round(e);
    }

    void merge(const HyperLogLog& other) {
        if (p != other.p) throw std::invalid_argument("Precision must match for merging");
        for (int i = 0; i < m; ++i) {
            if (registers[i] < other.registers[i]) {
                registers[i] = other.registers[i];
            }
        }
    }

private:

    int p;
    int m;
    std::vector<uint8_t> registers;
    double alpha;

    void compute_alpha() {
        if (m == 16) alpha = 0.673;
        else if (m == 32) alpha = 0.697;
        else if (m == 64) alpha = 0.709;
        else alpha = 0.7213 / (1 + 1.079 / m);
    }

    void process_hash(uint64_t hash) {
        uint32_t index = hash >> (64 - p);
        uint64_t w = hash & ((1ULL << (64 - p)) - 1);
        uint8_t rho_val = w == 0 ? (64 - p + 1) : (clz(w) - p + 1);
        if (registers[index] < rho_val) {
            registers[index] = rho_val;
        }
    }

};