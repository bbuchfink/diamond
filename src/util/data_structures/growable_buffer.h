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

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>

template<typename T = char>
class GrowableBuffer {
public:
    GrowableBuffer() = default;

    explicit GrowableBuffer(std::size_t initial_capacity) {
        reserve_discarding(initial_capacity);
    }

    GrowableBuffer(const GrowableBuffer&) = delete;
    GrowableBuffer& operator=(const GrowableBuffer&) = delete;

    GrowableBuffer(GrowableBuffer&&) noexcept = default;
    GrowableBuffer& operator=(GrowableBuffer&&) noexcept = default;

    T* data() noexcept {
        return buffer_.get();
    }

    const T* data() const noexcept {
        return buffer_.get();
    }

    std::size_t capacity() const noexcept {
        return capacity_;
    }

    bool empty() const noexcept {
        return capacity_ == 0;
    }

    void clear() noexcept {
    }

    void ensure_capacity(std::size_t min_capacity) {
        if (min_capacity <= capacity_) {
            return;
        }

        reserve_discarding(grown_capacity(min_capacity));
    }

    void reserve_discarding(std::size_t new_capacity) {
        if (new_capacity <= capacity_) {
            return;
        }
        std::unique_ptr<T[]> new_buffer(new T[new_capacity]);
        buffer_ = std::move(new_buffer);
        capacity_ = new_capacity;
    }

private:
    static std::size_t grown_capacity(std::size_t min_capacity) {
        constexpr std::size_t max = std::numeric_limits<std::size_t>::max();

        std::size_t new_capacity = 64;

        while (new_capacity < min_capacity) {
            if (new_capacity > max / 2) {
                if (min_capacity > max) {
                    throw std::length_error("buffer capacity overflow");
                }
                return min_capacity;
            }

            new_capacity *= 2;
        }

        return new_capacity;
    }

    std::unique_ptr<T[]> buffer_;
    std::size_t capacity_ = 0;
};