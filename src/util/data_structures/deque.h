#pragma once

#include <list>
#include <vector>
#include "../../basic/config.h"

template<typename _t>
struct Deque {

	typedef std::vector<_t> Bucket;

	Deque():
		max_size_(config.deque_bucket_size / sizeof(_t))
	{
		buckets.emplace_back();
		buckets.back().reserve(max_size_);
	}

	void push_back(const _t* ptr, size_t n) {
		if (buckets.back().size() + n > max_size_) {
			buckets.emplace_back();
			buckets.back().reserve(max_size_);
		}
		buckets.back().insert(buckets.back().end(), ptr, ptr + n);
	}

	size_t size() const {
		size_t n = 0;
		for (const Bucket& b : buckets)
			n += b.size();
		return n;
	}

	void move(std::vector<_t>& dst) {
		if (buckets.size() == 1 && dst.empty())
			dst = std::move(buckets.front());
		else {
			for (const Bucket& b : buckets)
				dst.insert(dst.end(), b.begin(), b.end());
		}
		buckets.clear();
	}

private:

	std::list<Bucket> buckets;
	const size_t max_size_;

};