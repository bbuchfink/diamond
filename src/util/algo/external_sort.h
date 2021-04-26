/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include <vector>
#include <numeric>
#include <queue>
#include "../io/temp_file.h"
#include "../io/input_file.h"
#define _REENTRANT
#include "../lib/ips4o/ips4o.hpp"
#include "../basic/config.h"

template<typename _t>
static inline size_t alloc_size(const _t& x) {
	return sizeof(_t);
}

static inline size_t alloc_size(const std::pair<std::string, uint32_t>& x) {
	return 4 + sizeof(std::string) + x.first.length();
}

template<typename Type, typename Cmp = std::less<Type>>
struct ExternalSorter {

	typedef Type Value;

	static const size_t BUCKET_SIZE = 2 * GIGABYTES;

	ExternalSorter(Cmp cmp = Cmp()) :
		cmp_(cmp),
		count_(0),
		size_(0)
	{
	}

	void push(const Type& x) {
		++count_;
		size_ += alloc_size(x);
		buf_.push_back(x);
		if (size_ > BUCKET_SIZE)
			flush();
	}

	void init_read() {
		flush();
		Entry e;
		for (size_t i = 0; i < files_.size(); ++i) {
			if(get_entry(i, e))
				queue_.push(e);
		}
	}

	bool good() const {
		return !queue_.empty();
	}

	const Type& operator*() const {
		return queue_.top().value;
	}

	void operator++() {
		const size_t b = queue_.top().bucket;
		queue_.pop();
		Entry e;
		if (get_entry(b, e))
			queue_.push(e);
	}

	size_t count() const {
		return count_;
	}

private:

	struct Entry {
		size_t bucket;
		Type value;
		bool operator<(const Entry& e) const {
			return Cmp()(e.value, value);
		}
	};

	struct CmpIdx {
		bool operator()(uint32_t x, uint32_t y) const {
			return Cmp()(data_[x], data_[y]);
		}
		typename std::vector<Type>::const_iterator data_;
	};

	bool get_entry(size_t bucket, Entry& e) {
		try {
			(*files_[bucket]) >> e.value;
		}
		catch (EndOfStream&) {
			files_[bucket]->close_and_delete();
			delete files_[bucket];
			return false;
		}
		e.bucket = bucket;
		return true;
	}
	
	void flush() {
		std::vector<uint32_t> idx_(buf_.size());
		std::iota(idx_.begin(), idx_.end(), 0);
		ips4o::parallel::sort(idx_.begin(), idx_.end(), CmpIdx{ buf_.begin() }, config.threads_);
		TempFile f;
		for (uint32_t i : idx_)
			f << buf_[i];
		files_.push_back(new InputFile(f));
		buf_.clear();
		size_ = 0;
	}

	const Cmp cmp_;
	size_t count_;
	size_t size_;
	std::vector<InputFile*> files_;
	std::vector<Type> buf_;
	std::priority_queue<Entry> queue_;

};