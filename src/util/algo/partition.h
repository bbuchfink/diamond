#pragma once
#include <stddef.h>

template<typename T = size_t>
struct Partition
{
	Partition() : parts(0), items(0), size_(0), remainder(0)
	{ }
	Partition(T items, T parts) : parts(std::min(parts, items)), items(items)
	{
		if (this->parts > 0) {
			size_ = items / this->parts;
			remainder = items % this->parts;
		}
		else {
			size_ = 0;
			remainder = 0;
		}
	}
	T begin(T i) const
	{
		T b = std::min(i, remainder); return b * (size_ + 1) + (i - b) * size_;
	}
	T end(T i) const
	{
		return begin(i) + size(i);
	}
	T size(T i) const
	{
		return i < remainder ? (size_ + 1) : size_;
	}
	const T parts;
private:
	T items, size_, remainder;
};
