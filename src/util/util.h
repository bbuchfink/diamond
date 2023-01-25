/****
DIAMOND protein aligner
Copyright (C) 2013-2022 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include <assert.h>
#include <vector>
#include <algorithm>
#include <string>
#include <limits>
#include <set>
#include <mutex>
#include <math.h>

template<typename T>
inline T div_up(T x, T m)
{
	return (x + (m - 1)) / m;
}

template<typename T>
inline T round_up(T x, T m)
{
	return div_up(x, m) * m;
}

std::vector<std::string> tokenize(const char* str, const char* delimiters);
extern const char dir_separator;
std::string extract_dir(const std::string &s);

struct Sd
{
	Sd():
		A(0),
		Q(0),
		k(1)
	{}
	Sd(const std::vector<Sd> &groups);
	void add(double x)
	{
		const double d = x - A;
		Q += (k - 1) / k*d*d;
		A += d / k;
		++k;
	}
	double mean() const
	{
		return A;
	}
	double sd() const
	{
		return sqrt(Q / (k - 1));
	}
private:
	double A, Q, k;	
};

std::string to_upper_case(const std::string& s);
std::string to_lower_case(const std::string& s);

template<typename T>
struct Matrix
{
	Matrix(int64_t rows, int64_t cols) :
		cols_(cols),
		data_(rows* cols, 0)
	{}
	T* operator[](int64_t i)
	{
		return &data_[i * cols_];
	}
private:
	int64_t cols_;
	std::vector<T> data_;
};

template<int n>
inline int round_down(int x)
{
	return (x / n) * n;
}

template<int n>
inline int round_up(int x)
{
	return ((x + n - 1) / n) * n;
}

std::string print_char(char c);

template<typename T1, typename T2>
T1 percentage(T2 x, T2 y)
{
	return x * (T1)100 / y;
}

void print_binary(uint64_t x);

template<typename T1, typename T2>
inline typename std::enable_if<std::is_signed<T2>::value, T1>::type safe_cast(T2 x)
{
	if (x > (T2)std::numeric_limits<T1>::max() || x < (T2)std::numeric_limits<T1>::min())
		throw std::runtime_error("Integer value out of bounds.");
	return (T1)x;
}

template<typename T1, typename T2>
inline typename std::enable_if<std::is_unsigned<T2>::value, T1>::type safe_cast(T2 x)
{
	if (x > (T2)std::numeric_limits<T1>::max())
		throw std::runtime_error("Integer value out of bounds.");
	return (T1)x;
}


struct IndexIterator
{
	IndexIterator(size_t i) :
		i(i)
	{}
	size_t operator*() const
	{
		return i;
	}
	bool operator!=(const IndexIterator &rhs) const
	{
		return i != rhs.i;
	}
	IndexIterator& operator++()
	{
		++i;
		return *this;
	}
	size_t i;
};

inline double megabytes(size_t x)
{
	return (double)x / (1 << 20);
}

template<typename T>
inline T make_multiple(T x, T m)
{
	if (x % m == 0)
		return x;
	T d = m - x % m;
	if (std::numeric_limits<T>::max() - d < x)
		return x;
	return x + d;
}

std::string hex_print(const char* x, int len);
std::set<int32_t> parse_csv(const std::string& s);
std::string join(const char *c, const std::vector<std::string> &v);

template<typename T, typename F>
auto apply(const std::vector<T> &v, F f) -> std::vector<typename std::result_of<F(T)>::type> {
	std::vector<typename std::result_of<F(T)>::type> r;
	r.reserve(v.size());
	for (const auto &i : v)
		r.push_back(f(i));
	return r;
}

template<typename T1, typename T2>
std::vector<std::pair<T1, T2>> combine(const std::vector<T1> &v1, const std::vector<T2> &v2) {
	std::vector<std::pair<T1, T2>> r;
	r.reserve(v1.size());
	for (size_t i = 0; i < v1.size(); ++i)
		r.emplace_back(v1[i], v2[i]);
	return r;
}

template<typename It, typename Key>
struct AsyncKeyMerger {
	AsyncKeyMerger(It begin, It end, const Key& key) :
		it_(begin),
		end_(end),
		key_(key)
	{}
	std::pair<It, It> operator++() {
		std::lock_guard<std::mutex> lock(mtx_);
		if (it_ == end_)
			return { end_, end_ };
		It begin = it_++;
		auto key = key_(*begin);
		while (it_ != end_ && key_(*it_) == key) ++it_;
		return { begin, it_ };
	}
private:
	It it_;
	const It end_;
	const Key key_;
	std::mutex mtx_;
};

template<typename It, typename Key>
struct KeyMergeIterator {
	using value_type = typename std::iterator_traits<It>::value_type;
	typedef typename std::result_of<Key(value_type&)>::type KeyType;
	KeyMergeIterator(const It& begin, const It& end, const Key& key) :
		end_(end),
		begin_(begin),
		key_end_(begin),
		get_key_(key)
	{
		if (begin == end)
			return;
		next_key_ = key(*begin);
		assert(begin != end);
		this->operator++();
	}
	void operator++()
	{
		begin_ = key_end_;
		if (begin_ == end_)
			return;
		key_ = next_key_;
		++key_end_;
		while (key_end_ != end_ && (next_key_ = get_key_(*key_end_)) == key_) ++key_end_;
	}
	bool good() const
	{
		return begin_ != end_;
	}
	It& begin()
	{
		return begin_;
	}
	It& end()
	{
		return key_end_;
	}
	KeyType key() const {
		return key_;
	}
	size_t count() const {
		return size_t(key_end_ - begin_);
	}
private:
	const It end_;
	It begin_, key_end_;
	const Key get_key_;
	KeyType key_, next_key_;
};

template<typename It, typename Key>
KeyMergeIterator<It, Key> inline merge_keys(const It& begin, const It& end, const Key& key) {
	return KeyMergeIterator<It, Key>(begin, end, key);
}