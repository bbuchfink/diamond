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
