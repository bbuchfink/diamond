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
#include <algorithm>
#include <tuple>
#include <vector>
#include <numeric>
#include <queue>
#include "../io/temp_file.h"
#include "../io/input_file.h"

Serializer& operator<<(Serializer& s, const std::tuple<std::string, uint32_t>& x);

template<typename _t, size_t _key>
struct ExternalSorter {

	typedef typename std::tuple_element<_key, _t>::type Key;

	static const size_t BUCKET_SIZE = 0x80000000;

	ExternalSorter() :
		bucket_size_(std::max(BUCKET_SIZE/sizeof(_t), (size_t)1)),
		count_(0)
	{
	}

	void push(const _t& x) {
		++count_;
		buf_.push_back(x);
		if (buf_.size() > bucket_size_)
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

	const _t& get() {
		return queue_.top().value;
	}

	void operator++() {
		const size_t b = queue_.top().bucket;
		queue_.pop();
		Entry e;
		if (get_entry(b, e))
			queue_.push(e);
	}

	static ExternalSorter end() {
		return ExternalSorter{ 0,0 };
	}

	bool operator==(const ExternalSorter& s) {
		return !good() && s.bucket_size_ == 0 ? true : false;
	}

private:

	struct Entry {
		size_t bucket;
		_t value;
		bool operator<(const Entry& e) const {
			return std::get<_key>(value) > std::get<_key>(e.value);
		}
	};

	struct Cmp {
		bool operator()(uint32_t x, uint32_t y) const {
			return std::get<_key>(data_[x]) < std::get<_key>(data_[y]);
		}
		typename std::vector<_t>::const_iterator data_;
	};

	bool get_entry(size_t bucket, Entry& e) {
		try {
			files_[bucket] >> e.value;
		}
		catch (EndOfStream&) {
			files_[bucket].close_and_delete();
			return false;
		}
		e.bucket = bucket;
		return true;
	}
	
	void flush() {
		std::vector<uint32_t> idx_(buf_.size());
		std::iota(idx_.begin(), idx_.end(), 0);
		std::sort(idx_.begin(), idx_.end(), Cmp{ buf_.begin() });
		TempFile f;
		for (uint32_t i : idx_)
			f << buf_[i];
		files_.emplace_back(f);
		buf_.clear();
	}

	const size_t bucket_size_;
	size_t count_;
	std::vector<InputFile> files_;
	std::vector<_t> buf_;
	std::priority_queue<Entry> queue_;

};

template<typename _t1, typename _t2, size_t _key1, size_t _key2>
struct TupleJoinIterator {

	typedef ExternalSorter<_t1, _key1> Sorter1;
	typedef ExternalSorter<_t2, _key2> Sorter2;
	typedef typename Sorter1::Key Key;

	TupleJoinIterator(Sorter1& sorter1, Sorter2& sorter2):
		sorter1_(sorter1),
		sorter2_(sorter2)
	{
		sorter1.init_read();
		sorter2.init_read();
	}
	
	bool good() {
		return sorter1_.good() && sorter2_.good();
	}

	std::pair<_t1, _t2> operator*() const {
		return { sorter1_.get(), sorter2_.get() };
	}

	void operator++() {
		++sorter1_;
		++sorter2_;
		if (!good())
			return;
		auto key1 = std::get<_key1>(sorter1_.get());
		auto key2 = std::get<_key2>(sorter2_.get());
		do {
			if (key1 < key2) {
				++sorter1_;
				key1 = std::get<_key1>(sorter1_.get());
			}
			else if (key2 < key1) {
				++sorter2_;
				key2 = std::get<_key2>(sorter2_.get());
			}
			else
				return;
		} while (good());
	}

private:

	Sorter1& sorter1_;
	Sorter2& sorter2_;

};