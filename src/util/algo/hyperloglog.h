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

        if (all_zero) return 0;

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
        uint32_t index = uint32_t(hash >> (64 - p));
        uint64_t w = hash & ((1ULL << (64 - p)) - 1);
        uint8_t rho_val = w == 0 ? (64 - p + 1) : (clz(w) - p + 1);
        if (registers[index] < rho_val) {
            registers[index] = rho_val;
        }
    }

};