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
#include <cstdint>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

class DegreePartition {
public:
    struct Bucket {
        uint64_t first_degree;
        uint64_t last_degree;
        uint64_t node_count;
    };

    DegreePartition(const std::unordered_map<uint64_t, uint64_t>& degrees, size_t desired_bucket_count) {
        if (desired_bucket_count == 0) {
            throw std::invalid_argument("desired_bucket_count must be > 0");
        }

        std::vector<std::pair<uint64_t, uint64_t>> hist;
        hist.reserve(degrees.size());

        for (const auto& i : degrees) {
            if (i.second != 0) {
                hist.emplace_back(i.first, i.second);
            }
        }

        std::sort(hist.begin(), hist.end(),
            [](const std::pair<uint64_t, uint64_t>& a, const std::pair<uint64_t, uint64_t>& b) {
                return a.first < b.first;
            });

        if (hist.empty()) {
            return;
        }

        const std::size_t n = hist.size();
        const std::size_t bucket_count = std::min(desired_bucket_count, n);

        std::vector<uint64_t> prefix(n + 1, 0);

        for (std::size_t i = 0; i < n; ++i) {
            prefix[i + 1] = checked_add(prefix[i], hist[i].second);
        }

        total_nodes_ = prefix.back();

        std::size_t start = 0;

        for (std::size_t b = 0; b < bucket_count; ++b) {
            std::size_t end_prefix_index;

            if (b + 1 == bucket_count) {
                end_prefix_index = n;
            }
            else {
                const std::size_t remaining_buckets = bucket_count - b - 1;

                const std::size_t min_p = start + 1;
                const std::size_t max_p = n - remaining_buckets;

                const long double ideal_cumulative =
                    static_cast<long double>(total_nodes_) *
                    static_cast<long double>(b + 1) /
                    static_cast<long double>(bucket_count);

                auto begin_it = prefix.begin() + static_cast<std::ptrdiff_t>(min_p);
                auto end_it = prefix.begin() + static_cast<std::ptrdiff_t>(max_p + 1);

                auto it = std::lower_bound(
                    begin_it,
                    end_it,
                    ideal_cumulative,
                    [](uint64_t lhs, long double rhs) {
                        return static_cast<long double>(lhs) < rhs;
                    }
                );

                if (it == end_it) {
                    end_prefix_index = max_p;
                }
                else {
                    end_prefix_index =
                        static_cast<std::size_t>(std::distance(prefix.begin(), it));
                }

                if (end_prefix_index > min_p) {
                    const std::size_t prev = end_prefix_index - 1;

                    if (distance(prefix[prev], ideal_cumulative) <=
                        distance(prefix[end_prefix_index], ideal_cumulative)) {
                        end_prefix_index = prev;
                    }
                }
            }

            const uint64_t first_degree = hist[start].first;

            const uint64_t last_degree =
                end_prefix_index < n
                ? hist[end_prefix_index].first - 1
                : hist.back().first;

            const uint64_t node_count =
                prefix[end_prefix_index] - prefix[start];

            buckets_.push_back(Bucket{
                first_degree,
                last_degree,
                node_count
                });

            start = end_prefix_index;
        }
    }

    const std::vector<Bucket>& buckets() const {
        return buckets_;
    }

    uint64_t total_nodes() const {
        return total_nodes_;
    }

    std::size_t size() const {
        return buckets_.size();
    }

    std::size_t bucket_index(uint64_t degree) const {
        auto it = std::upper_bound(
            buckets_.begin(),
            buckets_.end(),
            degree,
            [](uint64_t degree, const Bucket& bucket) {
                return degree < bucket.first_degree;
            }
        );

        if (it == buckets_.begin()) {
            throw std::out_of_range("degree is below partition range");
        }

        --it;

        if (degree > it->last_degree) {
            throw std::out_of_range("degree is above partition range");
        }

        return static_cast<std::size_t>(std::distance(buckets_.begin(), it));
    }

    const Bucket& bucket_for_degree(uint64_t degree) const {
        return buckets_.at(bucket_index(degree));
    }

private:
    static uint64_t checked_add(uint64_t a, uint64_t b) {
        if (b > std::numeric_limits<uint64_t>::max() - a) {
            throw std::overflow_error("node count overflow");
        }
        return a + b;
    }

    static long double distance(uint64_t value, long double target) {
        const long double diff = static_cast<long double>(value) - target;
        return diff < 0 ? -diff : diff;
    }

    std::vector<Bucket> buckets_;
    uint64_t total_nodes_ = 0;
};