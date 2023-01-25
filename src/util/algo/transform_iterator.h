#pragma once
#include <type_traits>

template<typename It, typename F>
struct TransformIterator {
	using difference_type = typename It::difference_type;
	using iterator_category = std::random_access_iterator_tag;
	using reference = typename std::result_of<F(const typename It::value_type&)>::type;
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