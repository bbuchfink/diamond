#ifndef DOUBLE_ITERATOR_H_
#define DOUBLE_ITERATOR_H_

#include <utility>

template<typename _it1, typename _it2, typename _t1, typename _t2>
struct DoubleIterator {
	typedef ptrdiff_t difference_type;
	std::pair<_t1, _t2> operator*() {
		return { *it1, *it2 };
	}
	difference_type operator-(const DoubleIterator &x) const {
		return it1 - x.it1;
	}
	DoubleIterator operator-(ptrdiff_t d) {
		return { it1 - d,it2 - d };
	}
	_it1 it1;
	_it2 it2;
};

#endif