#ifndef DOUBLE_ARRAY_H_
#define DOUBLE_ARRAY_H_

#include <stdint.h>
#include <stddef.h>
#include <algorithm>
#include "../range.h"

template<typename _t>
struct DoubleArray {

	DoubleArray()
	{}

	DoubleArray(void* data, size_t size = 0):
		data_((char*)data),
		size_(size)
	{}

	struct Iterator {

		Iterator(void *ptr) :
			ptr_((char*)ptr),
			end_(nullptr)
		{}

		Iterator(char *ptr, char *end) :
			ptr_(ptr),
			end_(end)
		{}

		uint32_t& count() {
			return *(uint32_t*)ptr_;
		}

		Range<_t*> operator*() {
			return { (_t*)(ptr_ + 4), (_t*)(ptr_ + 4) + count() };
		}

		Range<_t*>* operator->() {
			range_ = this->operator*();
			return &range_;
		}

		Iterator& operator++() {
			next();
			while (ptr_ < end_ && count() == 0)
				ptr_ += (*(uint32_t*)(ptr_ + 4)) * sizeof(_t) + 4;
			return *this;
		}

		void next() {
			ptr_ += count() * sizeof(_t) + 4;
		}

		operator bool() const {
			return ptr_ < end_;
		}

		void erase() {
			*(uint32_t*)(ptr_ + 4) = count();
			count() = 0;
		}

		ptrdiff_t operator-(const Iterator &x) const {
			return ptr_ - x.ptr_;
		}

		bool operator==(const Iterator &x) const {
			return ptr_ == x.ptr_;
		}

	private:

		Range<_t*> range_;
		char *ptr_, *end_;

		friend struct DoubleArray;

	};

	Iterator begin() {
		return Iterator(data_, data_ + size_);
	}

	void set_end(const Iterator &it) {
		size_ = it.ptr_ - data_;
	}

	void append(const DoubleArray &d) {
		if (d.data_ == data_ + size_) {
			size_ += d.size_;
			return;
		}
		std::copy(d.data_, d.data_ + d.size_, data_ + size_);
		size_ += d.size_;
	}

	uint32_t offset(const Iterator &it) const {
		return uint32_t(it.ptr_ - data_);
	}

	_t& operator[](uint32_t i) {
		return *(_t*)(data_ + i);
	}

private:

	char *data_;
	size_t size_;

};

#endif