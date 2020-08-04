/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
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
#include <vector>
#include <algorithm>
#include <iostream>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <limits>
#include <stdexcept>
#include <stdio.h>
#include <stdint.h>
#include <set>
#include <random>
#include "simd.h"
#include "../basic/const.h"

using std::vector;
using std::string;
using std::set;

template<typename _t=size_t>
struct partition
{
	_t items, parts, size, remainder;
	partition() : items(0), parts(0), size(0), remainder(0)
	{ }
	partition(_t items, _t parts) : items(items), parts(std::min(parts, items))
	{
		if(this->parts > 0) {
			size = items / this->parts;
			remainder = items % this->parts;
		} else {
			size = 0;
			remainder = 0;
		}
	}
	_t getMin(_t i) const
	{ _t b = std::min(i, remainder); return b*(size+1) + (i-b)*size; }
	_t getMax(_t i) const
	{ return getMin(i) + getCount(i); }
	_t getCount(_t i) const
	{ return i < remainder ? (size + 1) : size; }
};

#ifdef __SSE2__
inline void print(const __m128i &x)
{
	char *p=(char*)&x;
	for(unsigned i=0;i<16;++i)
		std::cout << int(*(p++)) << ' ';
	std::cout << std::endl;
}
#endif

template<typename _it, typename _key>
inline vector<size_t> map_partition(_it begin, _it end, const _key& key, size_t min_size, size_t max_segments, size_t min_segments)
{
	const size_t n = end - begin;
	const ::partition<size_t> p (n, std::max(min_segments, std::min(max_segments, n/min_size)));
	vector<size_t> v (p.parts+1);
	v[0] = p.getMin(0);
	v[p.parts] = p.getMax(p.parts-1);
	for(unsigned i=0;i<p.parts-1;++i) {
		size_t e = p.getMax(i);
		if(v[i] >= e) {
			v[i+1] = v[i];
			continue;
		}
		while(e < n && key(*(begin+e)) == key(*(begin+e-1)))
			++e;
		v[i+1] = e;
	}
	return v;
}

template<typename _t>
inline _t div_up(_t x, _t m)
{ return (x + (m-1)) / m; }

template<typename _t>
inline _t round_up(_t x, _t m)
{ return div_up(x, m) * m; }

inline vector<string> tokenize(const char *str, const char *delimiters)
{
	vector<string> out;
	string token;
	while(*str != 0) {
		while(*str != 0 && strchr(delimiters, *str))
			++str;
		token.clear();
		while(*str != 0 && strchr(delimiters, *str) == nullptr)
			token += *(str++);
		if(token.length() > 0)
			out.push_back(token);
	}
	if(out.size() == 0)
		out.push_back(string ());
	return out;
}

template<typename _t1, typename _t2>
struct Pair
{
	Pair():
		first (),
		second ()
	{ }
	Pair(const _t1 &first, const _t2 &second):
		first (first),
		second (second)
	{ }
	bool operator<(const Pair &rhs) const
	{
		return first < rhs.first;
	}
	_t1 first;
	_t2 second;
};

inline size_t find_first_of(const char *s, const char *delimiters)
{
	const char *t = s;
	while(*t && strchr(delimiters, *t) == 0)
		++t;
	return t-s;
}

inline string blast_id(const string &title)
{
	return title.substr(0, find_first_of(title.c_str(), Const::id_delimiters));
}

inline void get_title_def(const string &s, string &title, string &def)
{
	const size_t i = find_first_of(s.c_str(), Const::id_delimiters);
	title = s.substr(0, i);
	if (i >= s.length())
		def.clear();
	else
		def = s.substr(i + 1);
}

inline size_t print_str(char* buf, const char *s, size_t n)
{
	memcpy(buf, s, n);
	*(buf+n) = 0;
	return n;
}

inline size_t print_str(char *buf, const char *s, const char *delimiters)
{ return print_str(buf, s, find_first_of(s, delimiters)); }

inline string* get_str(const char *s, const char *delimiters)
{
	return new string (s, find_first_of(s, delimiters));
}

template<typename _t, unsigned d1, unsigned d2>
struct Static_matrix
{
	_t* operator[](size_t i)
	{ return data_[i]; }
private:
	_t data_[d1][d2];
};

extern const char dir_separator;
string extract_dir(const string &s);

inline std::ostream& indent(std::ostream &str, unsigned n)
{
	for (unsigned i = 0; i < n; ++i)
		str << ' ';
	return str;
}

inline int abs_diff(unsigned x, unsigned y)
{
	return abs((int)x - int(y));
}

template<typename _t, int _n>
struct Static_vector
{
	Static_vector():
		n(0)
	{}
	_t& operator[](int i)
	{
		return data[i];
	}
	const _t& operator[](int i) const
	{
		return data[i];
	}
	int size() const
	{
		return n;
	}
	void push_back(const _t &x)
	{
		data[n++] = x;
	}
	void erase(int i)
	{
		memmove(&data[i], &data[i + 1], (--n - i)*sizeof(_t));
	}
private:
	_t data[_n];
	int n;
};

template<typename _t>
inline void assign_ptr(_t& dst, _t *src)
{
	dst = *src;
	delete src;
}

struct Sd
{
	Sd():
		A(0),
		Q(0),
		k(1)
	{}
	Sd(const vector<Sd> &groups);
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

inline string to_upper_case(const string &s)
{
	string r;
	for (string::const_iterator i = s.begin(); i != s.end(); ++i)
		r.push_back(toupper(*i));
	return r;
}

inline string to_lower_case(const string &s)
{
	string r;
	for (string::const_iterator i = s.begin(); i != s.end(); ++i)
		r.push_back(tolower(*i));
	return r;
}

template<typename _t>
struct Matrix
{
	void init(int rows, int cols, bool reset = false)
	{
		cols_ = cols;
		if (!reset) {
			data_.clear();
			data_.resize(rows*cols);
		}
		else {
			data_.clear();
			data_.insert(data_.begin(), rows*cols, _t());
		}
	}
	_t* operator[](int i)
	{
		return &data_[i*cols_];
	}
private:
	int cols_;
	vector<_t> data_;
};

template<int n>
inline int round_down(int x)
{
	return (x / n)*n;
}

template<int n>
inline int round_up(int x)
{
	return ((x + n - 1) / n)*n;
}

template<typename _t>
bool equal(const _t *ptr, unsigned n)
{
	const _t v = *ptr;
	const _t* end = (ptr++) + n;
	for (; ptr < end; ++ptr)
		if (*ptr != v)
			return false;
	return true;
}

inline string print_char(char c)
{
	char buf[16];
	if (c < 32)
		sprintf(buf, "ASCII %u", (unsigned)c);
	else
		sprintf(buf, "%c", c);
	return string(buf);
}

template<typename _t, int n>
struct Top_list
{
	_t& add(const _t &x)
	{
		for (int i = 0; i < n; ++i)
			if ((int)x >(int)data_[i]) {
				if (i < n - 1)
					memmove(&data_[i + 1], &data_[i], sizeof(data_)/n*(n - 1 - i));
				data_[i] = x;
				return data_[i];
			}
	}
	const _t& operator[](unsigned i) const
	{
		return data_[i];
	}
	_t& operator[](unsigned i)
	{
		return data_[i];
	}
	void sort()
	{
		std::sort(&data_[0], &data_[n]);
	}
private:
	_t data_[n];
};

template<typename _t1, typename _t2>
_t1 percentage(_t2 x, _t2 y)
{
	return x * (_t1)100 / y;
}

template<typename _t>
struct Numeric_vector : public std::vector<_t>
{
	Numeric_vector(size_t n):
		vector<_t>(n)
	{}
	Numeric_vector& operator+=(Numeric_vector &x)
	{
		for (size_t i = 0; i < this->size(); ++i)
			this->operator[](i) += x[i];
		return *this;
	}
	Numeric_vector& operator/=(double x)
	{
		for (size_t i = 0; i < this->size(); ++i)
			this->operator[](i) /= x;
		return *this;
	}
	friend std::ostream& operator<<(std::ostream &s, const Numeric_vector &x)
	{
		for (size_t i = 0; i < x.size(); ++i)
			s << x[i] << std::endl;
		return s;
	}
};

template<unsigned n>
inline unsigned get_distribution(const double *p, std::minstd_rand0 &random_engine)
{
	const double x = std::uniform_real_distribution<double>(0.0, 1.0)(random_engine);
	double s = 0;
	for (unsigned i = 0; i < n; ++i) {
		s += p[i];
		if (x < s)
			return i;
	}
	return n - 1;
}

inline void print_binary(uint64_t x)
{
	for (unsigned i = 0; i < 64; ++i) {
		std::cout << (x & 1);
		x >>= 1;
	}
}

template<typename _t>
inline _t safe_cast(size_t x)
{
	if (x > (size_t)std::numeric_limits<_t>::max())
		throw std::runtime_error("Integer value out of bounds.");
	return (_t)x;
}

struct Index_iterator
{
	Index_iterator(size_t i) :
		i(i)
	{}
	size_t operator*() const
	{
		return i;
	}
	bool operator!=(const Index_iterator &rhs) const
	{
		return i != rhs.i;
	}
	Index_iterator& operator++()
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

inline int make_multiple(int x, int m)
{
	return x % m == 0 ? x : x + m - x%m;
}

inline string hex_print(const char *x, int len) {
	string out;
	char d[3];
	for (int i = 0; i < len; i++) {
		sprintf(d, "%02x", (unsigned char)x[i]);
		out += d;
	}
	return out;
}

inline std::set<unsigned> parse_csv(const std::string &s)
{
	std::set<unsigned> r;
	std::vector<std::string> t(tokenize(s.c_str(), ","));
	for (std::vector<std::string>::const_iterator i = t.begin(); i != t.end(); ++i)
		if(!i->empty()) r.insert(atoi(i->c_str()));
	return r;
}

std::string join(const char *c, const std::vector<std::string> &v);

template<typename _t, typename _f>
auto apply(const std::vector<_t> &v, _f f) -> std::vector<typename std::result_of<_f(_t)>::type> {
	std::vector<typename std::result_of<_f(_t)>::type> r;
	r.reserve(v.size());
	for (const auto &i : v)
		r.push_back(f(i));
	return r;
}

template<typename _t1, typename _t2>
std::vector<std::pair<_t1, _t2>> combine(const std::vector<_t1> &v1, const std::vector<_t2> &v2) {
	std::vector<std::pair<_t1, _t2>> r;
	r.reserve(v1.size());
	for (size_t i = 0; i < v1.size(); ++i)
		r.emplace_back(v1[i], v2[i]);
	return r;
}