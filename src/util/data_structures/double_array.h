#ifndef DOUBLE_ARRAY_H_
#define DOUBLE_ARRAY_H_

#include "../range.h"

template<typename _t>
struct DoubleArray {

	DoubleArray()
	{}

	DoubleArray(_t* data, size_t size):
		data_(data),
		size_(size)
	{}

	struct Iterator {

		Iterator(_t *ptr, _t *end) :
			ptr_(ptr),
			end_(end)
		{}

		Range<_t*> operator*() {
			return { ptr_ + 1, ptr_ + 1 + *ptr_ };
		}

		Range<_t*>* operator->() {
			range_ = this->operator*();
			return &range_;
		}

		Iterator& operator++() {
			ptr_ += *ptr_ + 1;
			while (ptr_ < end_ && *ptr_ == 0)
				ptr_ += ptr_[1] + 1;
			range_ = { ptr_ + 1, ptr_ + 1 + *ptr_ };
			return *this;
		}

		operator bool() const {
			return ptr_ < end_;
		}

		void erase() {
			ptr_[1] = *ptr_;
			*ptr_ = 0;
		}

	private:

		Range<_t*> range_;
		_t *ptr_, *end_;

	};

	Iterator begin() {
		return Iterator(data_, data_ + size_);
	}

private:

	_t *data_;
	size_t size_;

};

#endif