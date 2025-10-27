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