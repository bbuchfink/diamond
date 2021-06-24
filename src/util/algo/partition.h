#pragma once
#include <algorithm>

template<typename _t = size_t>
struct Partition
{
	Partition() : parts(0), items(0), size_(0), remainder(0)
	{ }
	Partition(_t items, _t parts) : parts(std::min(parts, items)), items(items)
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
	_t begin(_t i) const
	{
		_t b = std::min(i, remainder); return b * (size_ + 1) + (i - b) * size_;
	}
	_t end(_t i) const
	{
		return begin(i) + size(i);
	}
	_t size(_t i) const
	{
		return i < remainder ? (size_ + 1) : size_;
	}
	const _t parts;
private:
	_t items, size_, remainder;
};
