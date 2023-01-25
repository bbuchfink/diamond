#pragma once
#include <stdint.h>
#include <algorithm>

template<typename T>
struct Array {

	using Size = int64_t;

	Array():
		ptr_(nullptr),
		size_(0)
	{}

	Array(Size alloc_size) :
		ptr_(new T[alloc_size]),
		size_(0)
	{}

	Array& operator=(Array&& a) noexcept {
		ptr_ = a.ptr_;
		size_ = a.size_;
		a.ptr_ = nullptr;
		return *this;
	}

	~Array() {
		delete[] ptr_;
	}

	T* data() {
		return ptr_;
	}

	T* begin() {
		return ptr_;
	}

	T* end() {
		return ptr_ + size_;
	}

	Size size() const {
		return size_;
	}

	void assign(const T& x) {
		*ptr_ = x;
		size_ = 1;
	}

	template<typename It>
	void assign(It begin, It end) {
		std::copy(begin, end, ptr_);
		size_ = end - begin;
	}

	template<typename It>
	void assign_reversed(It begin, It end) {
		T* p = ptr_;
		for (It i = end - 1; i >= begin; --i)
			*p++ = *i;
		size_ = end - begin;
	}

	template<typename It>
	void push_back(It begin, It end) {
		std::copy(begin, end, this->end());
		size_ += end - begin;
	}

	template<typename It>
	void push_back_reversed(It begin, It end) {
		T* p = this->end();
		for (It i = end - 1; i >= begin; --i)
			*p++ = *i;
		size_ += end - begin;
	}

	void push_back(Size n, const T& value) {
		T* p = end();
		std::fill(p, p + n, value);
		size_ += n;
	}

private:
	T* ptr_;
	Size size_;

};