/****
Copyright © 2013-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-2-Clause

#pragma once
#include <cmath>
#include <assert.h>
#include <vector>
#include <string>
#include <limits>
#include <mutex>

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

template <class To, class From,
	typename std::enable_if<std::is_integral<From>::value && std::is_integral<To>::value, int>::type = 0>
inline To safe_cast(From value) {
	static_assert(std::is_integral<To>::value, "Destination must be an integral type");
	using ToT = To;
	using FromT = From;
#if __cplusplus >= 201703L
	if constexpr (std::is_signed<FromT>::value && std::is_signed<ToT>::value) {
#else
	if (std::is_signed<FromT>::value && std::is_signed<ToT>::value) {
#endif
		if (value < std::numeric_limits<ToT>::min() ||
			value > std::numeric_limits<ToT>::max()) {
			throw std::runtime_error("safe_cast: out of range (signed -> signed)");
		}
	}
#if __cplusplus >= 201703L
	else if constexpr (std::is_signed<FromT>::value && !std::is_signed<ToT>::value) {
#else
	else if (std::is_signed<FromT>::value && !std::is_signed<ToT>::value) {
#endif
		if (value < 0) {
			throw std::runtime_error("safe_cast: negative value to unsigned");
		}
		using UFrom = typename std::make_unsigned<FromT>::type;
		if (static_cast<UFrom>(value) > std::numeric_limits<ToT>::max()) {
			throw std::runtime_error("safe_cast: overflow (signed -> unsigned)");
		}
	}
#if __cplusplus >= 201703L
	else if constexpr (!std::is_signed<FromT>::value && std::is_signed<ToT>::value) {
#else
	else if (!std::is_signed<FromT>::value && std::is_signed<ToT>::value) {
#endif
		using UTo = typename std::make_unsigned<ToT>::type;
		if (value > static_cast<UTo>(std::numeric_limits<ToT>::max())) {
			throw std::runtime_error("safe_cast: overflow (unsigned -> signed)");
		}
	}
	else {
		if (value > std::numeric_limits<ToT>::max()) {
			throw std::runtime_error("safe_cast: overflow (unsigned -> unsigned)");
		}
	}
	return static_cast<ToT>(value);
}

template <class To, class From,
	typename std::enable_if<std::is_floating_point<From>::value && std::is_integral<To>::value, int>::type = 0>
inline To safe_cast(From value) {
	static_assert(std::is_integral<To>::value, "Destination must be an integral type");
	if (!std::isfinite(value)) {
		throw std::runtime_error("safe_cast: non-finite value (NaN/Inf)");
	}
#if __cplusplus >= 201703L
	if constexpr (std::is_unsigned<To>::value) {
#else
	if (std::is_unsigned<To>::value) {
#endif
		if (value < 0) {
			throw std::runtime_error("safe_cast: negative value to unsigned");
		}
	}
	using LD = long double;
	const LD lv = static_cast<LD>(value);
	const LD tmin = static_cast<LD>(std::numeric_limits<To>::min());
	const LD tmax = static_cast<LD>(std::numeric_limits<To>::max());
	if (lv < tmin || lv > tmax) {
		throw std::runtime_error("safe_cast: out of range (float -> int)");
	}
	return static_cast<To>(lv);
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
#if defined(__cpp_lib_is_invocable) && __cpp_lib_is_invocable >= 201703L
	using KeyType = typename std::invoke_result<Key, value_type&>::type;
#else
	using KeyType = typename std::result_of<Key(value_type&)>::type;
#endif
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