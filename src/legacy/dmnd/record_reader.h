/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <stdint.h>
#include <string.h>

struct Finish {};

struct DynamicRecordReader {

	DynamicRecordReader(File& d) :
		d_(d)
	{
		d.read(size_);
		size_ = big_endian_byteswap(size_);
	}

	DynamicRecordReader& operator>>(unsigned long long& x)
	{
		if (size_ >= sizeof(unsigned long long)) {
			d_.read(x);
			x = big_endian_byteswap(x);
			size_ -= sizeof(unsigned long long);
		}
		else
			x = 0;
		return *this;
	}

	DynamicRecordReader& operator>>(unsigned long& x)
	{
		if (size_ >= sizeof(unsigned long)) {
			d_.read(x);
			x = big_endian_byteswap(x);
			size_ -= sizeof(unsigned long);
		}
		else
			x = 0;
		return *this;
	}

	template<typename T>
	DynamicRecordReader& read(const T* ptr, size_t count)
	{
		const size_t s = count * sizeof(T);
		if (size_ >= s) {
			d_.read((void*)ptr, count * sizeof(T));
			size_ -= s;
		}
		else
			memset((void*)ptr, 0, s);
		return *this;
	}

	DynamicRecordReader& operator>>(int& x)
	{
		if (size_ >= sizeof(int)) {
			d_.read(x);
			x = big_endian_byteswap(x);
			size_ -= sizeof(int);
		}
		else
			x = 0;

		return *this;
	}

	void operator>>(const Finish&) {
		if (size_ == 0)
			return;
		char* buf = new char[size_];
		d_.read(buf, size_);
		delete[] buf;
		size_ = 0;
	}

private:
	File& d_;
	uint64_t size_;

};