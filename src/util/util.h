/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef UTIL_H_
#define UTIL_H_

#ifdef __MMX__
#include <mmintrin.h>
#endif

#ifdef __SSE__
#include <xmmintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <string.h>
#include <sys/stat.h>
#include "../basic/const.h"
#include "../basic/seed.h"
#include "../basic/value.h"
#include "../basic/packed_loc.h"

using std::vector;

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

inline Letter set_critical(Letter v)
{
	return static_cast<unsigned>(v) | 0x80;
}

inline Letter mask_critical(Letter v)
{
	return static_cast<unsigned>(v) & 0x7F;
}

inline bool get_critical(Letter v)
{
	return (static_cast<unsigned>(v) & 0x80) != 0;
}

inline unsigned filter_treshold(unsigned n)
{
	return config.hit_cap * 256 / n;
}

inline bool position_filter(Loc l, unsigned treshold, Packed_seed s)
{
	return ((l ^ s) & 0xff) < treshold;
}

struct interval
{
	interval():
		begin_ (0),
		end_ (0)
	{ }
	interval(size_t begin, size_t end):
		begin_ (begin),
		end_ (end)
	{ }
	size_t length() const
	{ return end_ > begin_ ? end_ - begin_ : 0; }
	size_t overlap(const interval &rhs) const
	{ return intersect(*this, rhs).length(); }
	bool includes(size_t p) const
	{ return p >= begin_ && p < end_; }
	friend inline interval intersect(const interval &lhs, const interval &rhs)
	{ return interval (std::max(lhs.begin_, rhs.begin_), std::min(lhs.end_, rhs.end_)); }
	friend std::ostream& operator<<(std::ostream &os, const interval &x)
	{ os << "[" << x.begin_ << ";" << x.end_ << "]"; return os; }
	size_t begin_, end_;
};

inline void print(const __m128i &x)
{
	char *p=(char*)&x;
	for(unsigned i=0;i<16;++i)
		std::cout << int(*(p++)) << ' ';
	std::cout << std::endl;
}

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

template<typename _t>
struct Tls
{
	Tls(_t *&ptr):
		ptr_ (init(ptr))
	{ }
	_t& operator*() const
	{ return *ptr_; }
	_t* operator->() const
	{ return ptr_; }
private:
	static _t* init(_t *&ptr)
	{
		if (ptr == 0)
			ptr = new _t;
		return ptr;
	}
	_t *const ptr_;
};

template<typename _t>
_t& get_tls(_t *&ptr)
{
	if (ptr == 0)
		ptr = new _t;
	return *ptr;
}

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

inline __m128i _mm_set(int a)
{
	int y = a << 8 | a;
	__m128i z = __m128i();
	z = _mm_insert_epi16 (z, y, 0);
	z = _mm_insert_epi16 (z, y, 1);
	z = _mm_insert_epi16 (z, y, 2);
	z = _mm_insert_epi16 (z, y, 3);
	z = _mm_insert_epi16 (z, y, 4);
	z = _mm_insert_epi16 (z, y, 5);
	z = _mm_insert_epi16 (z, y, 6);
	z = _mm_insert_epi16 (z, y, 7);
	return z;
}

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

inline string to_string(unsigned val)
{
	std::stringstream ss;
	ss << val;
	return ss.str();
}

string extract_dir(const string &s);

#endif /* UTIL_H_ */
