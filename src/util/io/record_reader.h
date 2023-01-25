#ifndef RECORD_READER_H_
#define RECORD_READER_H_

#include <stdint.h>
#include <string.h>
#include <../basic/value.h>
#include "deserializer.h"

struct Finish {};

struct DynamicRecordReader {

	DynamicRecordReader(Deserializer &d) :
		d_(d)
	{
		d.varint = false;
		d >> size_;
	}

	DynamicRecordReader& operator>>(unsigned long long &x)
	{
		if (size_ >= sizeof(unsigned long long)) {
			d_ >> x;
			size_ -= sizeof(unsigned long long);
		}
		else
			x = 0;
		return *this;
	}

	DynamicRecordReader& operator>>(unsigned long &x)
	{
		if (size_ >= sizeof(unsigned long)) {
			d_ >> x;
			size_ -= sizeof(unsigned long);
		}
		else
			x = 0;
		return *this;
	}

	template<typename _t>
	DynamicRecordReader& read(const _t *ptr, size_t count)
	{
		const size_t s = count * sizeof(_t);
		if (size_ >= s) {
			d_.read(ptr, count);
			size_ -= s;
		}
		else
			memset((void*)ptr, 0, s);
		return *this;
	}

	DynamicRecordReader& operator>>(int& x)
    {
        if (size_ >= sizeof(int)) {
            d_ >> x;
            size_ -= sizeof(int);
        }
        else
            x = 0;

		return *this;
	}

    void operator>>(const Finish&) {
        if (size_ == 0)
            return;
        char *buf = new char[size_];
        d_.read(buf, size_);
        delete[] buf;
        size_ = 0;
    }

private:
	Deserializer &d_;
	uint64_t size_;

};

#endif