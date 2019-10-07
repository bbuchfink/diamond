#ifndef ARRAY_H_
#define ARRAY_H_

#include <algorithm>

template<typename _t, size_t n>
struct Array
{
	Array()
	{
		std::fill(data_, data_ + n, _t());
	}
	_t& operator[](size_t i)
	{
		return data_[i];
	}
	const _t& operator[](size_t i) const
	{
		return data_[i];
	}
	_t* begin()
	{
		return data_;
	}
	const _t* begin() const
	{
		return data_;
	}
	_t* end()
	{
		return data_ + n;
	}
	_t data_[n];
};

#endif