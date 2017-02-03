/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <limits>
#include <stdexcept>
#include "simd.h"
#include "../basic/const.h"
#include "../basic/packed_loc.h"

using std::vector;
using std::string;

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

struct interval
{
	interval():
		begin_ (0),
		end_ (0)
	{ }
	interval(unsigned begin, unsigned end):
		begin_ (begin),
		end_ (end)
	{ }
	unsigned length() const
	{ return end_ > begin_ ? end_ - begin_ : 0; }
	unsigned overlap(const interval &rhs) const
	{ return intersect(*this, rhs).length(); }
	double overlap_factor(const interval &rhs) const
	{
		return (double)overlap(rhs) / (double)length();
	}
	bool includes(unsigned p) const
	{ return p >= begin_ && p < end_; }
	friend inline interval intersect(const interval &lhs, const interval &rhs)
	{ return interval (std::max(lhs.begin_, rhs.begin_), std::min(lhs.end_, rhs.end_)); }
	friend std::ostream& operator<<(std::ostream &os, const interval &x)
	{ os << "[" << x.begin_ << ";" << x.end_ << "]"; return os; }
	bool operator<(const interval &rhs) const
	{
		return begin_ < rhs.begin_;
	}
	unsigned begin_, end_;
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

template<typename _val, typename _dir>
inline void print_seq(const _val* s, const _dir& d)
{
	unsigned i=0;
	while(get_dir(s,i,d) != 0xff) {
		std::cout << mask_critical(get_dir(s,i,d));
		++i;
	}
}

#ifdef _MSC_VER
#define TLS_PTR __declspec(thread)
#else
#define TLS_PTR __thread
#endif

struct Ptr_wrapper_base
{
	virtual ~Ptr_wrapper_base()
	{}
};

template<typename _t>
struct Ptr_wrapper : public Ptr_wrapper_base
{
	Ptr_wrapper(_t *ptr) :
		ptr(ptr)
	{}
	virtual ~Ptr_wrapper()
	{
		delete ptr;
	}
	_t *ptr;
};

struct TLS
{
	template<typename _t>
	static _t& get(_t *&ptr)
	{
		if (ptr == 0) {
			ptr = new _t;
			if (ptr_ == 0)
				ptr_ = new vector<Ptr_wrapper_base*>;
			ptr_->push_back(new Ptr_wrapper<_t>(ptr));
		}
		return *ptr;
	}
	static void clear()
	{
		if (ptr_ == 0)
			return;
		for (vector<Ptr_wrapper_base*>::iterator i = ptr_->begin(); i != ptr_->end(); ++i)
			delete *i;
		delete ptr_;
	}
private:
	static TLS_PTR vector<Ptr_wrapper_base*> *ptr_;
};

inline vector<string> tokenize(const char *str, const char *delimiters)
{
	vector<string> out;
	while(*str != 0) {
		while(*str != 0 && strchr(delimiters, *str))
			++str;
		string token;
		while(*str != 0 && strchr(delimiters, *str) == 0)
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

#ifdef __SSE2__
inline __m128i _mm_set(char a)
{
	const int x = (int)a;
	return _mm_set1_epi32(x << 24 | x << 16 | x << 8 | x);
}
#endif

template<typename _t, unsigned d1, unsigned d2>
struct Static_matrix
{
	_t* operator[](size_t i)
	{ return data_[i]; }
private:
	_t data_[d1][d2];
};

inline void auto_append_extension(string &str, const char *ext)
{
	size_t l = strlen(ext);
	if(str.length() < l || (str.length() >= l && str.substr(str.length()-l, string::npos) != ext))
		str += ext;
}

inline bool check_dir(const string &path)
{
#ifndef _MSC_VER
	struct stat sb;
	return stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode);
#else
	return true;
#endif
}

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

template<typename _t, size_t n>
struct Array
{
	Array()
	{}
	_t& operator[](size_t i)
	{
		return data_[i];
	}
	const _t& operator[](size_t i) const
	{
		return data_[i];
	}
	_t* begin()
	{
		return data_;
	}
	const _t* begin() const
	{
		return data_;
	}
	_t* end()
	{
		return data_ + n;
	}
	_t data_[n];
};

inline int ctz(uint64_t x)
{
#ifdef _MSC_VER
	if (x)
		return (int)__popcnt64((x ^ (x - 1)) >> 1);
	else
		return CHAR_BIT * sizeof(x);
#else
	return __builtin_ctzll(x);
#endif
}

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
		return data_[i*cols_];
	}
private:
	int cols_;
	vector<_t> data_;
};

inline int get_idx(const char **a, size_t n, const char *s)
{
	for (size_t i = 0; i < n; ++i)
		if (strcmp(a[i], s) == 0)
			return (int)i;
	return -1;
}

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

template<typename _t>
inline string to_string(_t val)
{
	std::stringstream ss;
	ss << val;
	return ss.str();
}

inline string print_char(char c)
{
	std::stringstream ss;
	if (c < 32)
		ss << "ASCII " << (unsigned)c;
	else
		ss << c;
	return ss.str();
}

template<typename _t, int n>
struct Top_list
{
	void add(const _t &x)
	{
		for (int i = 0; i < n; ++i)
			if ((int)x >(int)data_[i]) {
				if (i < n - 1)
					memmove(&data_[i + 1], &data_[i], sizeof(data_)/n*(n - 1 - i));
				data_[i] = x;
				return;
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
inline unsigned get_distribution(const double *p)
{
	const double x = (double)rand() / RAND_MAX;
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

#endif /* UTIL_H_ */
