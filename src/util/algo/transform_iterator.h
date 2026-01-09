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
#include <type_traits>

template<typename It, typename F>
struct TransformIterator {
	using difference_type = typename It::difference_type;
	using iterator_category = std::random_access_iterator_tag;
#if __cplusplus >= 201703L
	using reference = typename std::invoke_result<F, const typename It::value_type&>::type;
#else
	using reference = typename std::result_of<F(const typename It::value_type&)>::type;
#endif
	using value_type = reference;
	using pointer = reference*;
	TransformIterator(It it, const F& f) :
		it_(it),
		f_(f)
	{}
	reference operator*() const {
		return f_(*it_);
	}
	TransformIterator& operator++() {
		++it_;
		return *this;
	}
	TransformIterator& operator--() {
		--it_;
		return *this;
	}
	TransformIterator& operator+=(difference_type d) {
		it_ += d;
		return *this;
	}
	bool operator!=(const TransformIterator& it) const {
		return it_ != it.it_;
	}
	TransformIterator operator+(difference_type n) const {
		return TransformIterator(it_ + n, f_);
	}
	difference_type operator+(const TransformIterator& it) const {
		return it_ + it.it_;
	}
	difference_type operator-(const TransformIterator& it) const {
		return it_ - it.it_;
	}
private:
	It it_;
	F f_;
};

template<typename It, typename F>
TransformIterator<It, F> transform(It it, const F& f) {
	return TransformIterator<It, F>(it, f);
}

template<typename It, typename F>
TransformIterator<It, F> operator+(typename TransformIterator<It, F>::difference_type n, const TransformIterator<It, F>& it) {
	return it + n;
}
