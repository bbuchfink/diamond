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

#include <type_traits>
#include <iterator>
#include <functional>
#include <cstddef>
#include <utility>

template <class RAIt>
using iter_cat_t = typename std::iterator_traits<RAIt>::iterator_category;

template <class RAIt, class Compare = std::less<typename std::remove_reference<decltype(*std::declval<RAIt>())>::type>>
void insertion_sort(RAIt first, RAIt last, Compare comp = Compare()) {
    static_assert(std::is_base_of<std::random_access_iterator_tag, iter_cat_t<RAIt>>::value,
                  "Requires RandomAccessIterators");

    using T = typename std::remove_reference<decltype(*first)>::type;
    
    const std::ptrdiff_t n = last - first;
    if (n < 2) return;

    RAIt min_it = first;
    for (RAIt it = first + 1; it < last; ++it) {
        if (comp(*it, *min_it)) min_it = it;
    }
    if (min_it != first) {
        T tmp = std::move(*first);
        *first = std::move(*min_it);
        *min_it = std::move(tmp);
    }

    for (RAIt i = first + 1; i < last; ++i) {
        T val = std::move(*i);
        RAIt j = i;
        while (comp(val, *(j - 1))) {
            *j = std::move(*(j - 1));
            --j;
        }
        *j = std::move(val);
    }
}