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
#include <set>

template<typename It1, typename It2, typename Cmp, typename Value>
struct SortedListJoiner {

	SortedListJoiner(It1& it1, It2& it2, Cmp cmp, Value value) :
		it1_(&it1),
		it2_(&it2),
		cmp_(cmp),
		value_(value)
	{
	}

	bool good() {
		return it1_->good() && it2_->good();
	}

	typename std::result_of<Value(const typename It1::Value&, const typename It2::Value&)>::type operator*() const {
		return value_(**it1_, **it2_);
	}

	void operator++() {
		const auto v1 = **it1_;
		const auto v2 = **it2_;
		++(*it1_);
		if (!it1_->good())
			return;

		if (!cmp_(v2, **it1_) && !cmp_(**it1_, v2))
			return;

		++(*it2_);
		if (!it2_->good())
			return;

		if (!cmp_(v1, **it2_) && !cmp_(**it2_, v1))
			throw std::runtime_error("Duplicate keys: " + v1.first);

		do {
			if (cmp_(**it1_, **it2_))
				++(*it1_);
			else if (cmp_(**it2_, **it1_))
				++(*it2_);
			else
				return;
		} while (good());
	}

private:

	It1* it1_;
	It2* it2_;
	const Cmp cmp_;
	const Value value_;

};

template<typename It1, typename It2, typename Cmp, typename Value>
SortedListJoiner<It1, It2, Cmp, Value> join_sorted_lists(It1& it1, It2& it2, Cmp cmp, Value value) {
	return SortedListJoiner<It1, It2, Cmp, Value>(it1, it2, cmp, value);
}

template<typename It, typename Key, typename Value>
struct KeyMerger
{

	typedef typename std::result_of<Key(const typename It::Value&)>::type KeyType;
	typedef typename std::result_of<Value(const typename It::Value&)>::type ValueType;

	KeyMerger(It& it, KeyType begin):
		it_(&it),
		key_(begin)
	{}

	KeyMerger& operator++() {
		++key_;
		return *this;
	}

	std::set<ValueType> operator*() {
		std::set<ValueType> v;
		while (it_->good() && Key()(**it_) == key_) {
			v.insert(Value()(**it_));
			++(*it_);
		}
		return v;
	}

	KeyType key() const {
		return key_;
	}

private:

	It* it_;
	KeyType key_;

};

template<typename It, typename Key, typename Value>
KeyMerger<It, Key, Value> merge_keys(It& it, Key key, Value value, typename std::result_of<Key(const typename It::Value&)>::type begin) {
	return KeyMerger<It, Key, Value>(it, begin);
}

template<typename Type1, typename Type2>
struct First {
	Type1 operator()(const std::pair<Type1, Type2>& p) const {
		return p.first;
	}
};

template<typename Type1, typename Type2>
struct Second {
	Type2 operator()(const std::pair<Type1, Type2>& p) const {
		return p.second;
	}
};