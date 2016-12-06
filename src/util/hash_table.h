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

#ifndef HASH_TABLE_H_
#define HASH_TABLE_H_

#include <exception>
#include <stdint.h>
#include <memory>
#include <algorithm>
#include <string.h>
#include "hash_function.h"

using std::auto_ptr;

struct hash_table_overflow_exception : public std::exception
{
	virtual const char* what() const throw()
	{
		return "Hash table overflow";
	}
};

template<typename T, T value> struct value_compare
{
bool operator()(T x) const
{
	return x == value;
}
};

template<typename _K, typename _V, typename _E, typename _H> class hash_table
{

public:

	struct entry
	{
		_K	key;
		_V	value;
	}; // __attribute__((packed));

	hash_table(size_t size) :
		table(new entry[size]),
		size_(size)
	{
		memset(table, 0, size_ * sizeof(entry));
	}

	~hash_table()
	{
		delete[] table;
	}

	entry* operator[](_K key) const
	{
		entry *entry = get_entry(key);
		if (_E()(entry->value))
			return NULL;
		return entry;
	}

	void insert(_K key, _V value)
	{
		entry *entry = get_entry(key);
		if (_E()(entry->value))
			entry->key = key;
		entry->value = value;
	}

	size_t size() const
	{
		return size_;
	}

	size_t count() const
	{
		size_t n(0);
		for (size_t i = 0; i < size_; ++i)
			if (!_E()(table[i].value))
				++n;
		return n;
	}

private:

	entry* get_entry(_K key) const
	{
		entry *p = &table[_H()(key) % size_];
		bool wrapped = false;
		while (p->key != key && !_E()(p->value)) {
			++p;
			if (p == &table[size_]) {
				if (wrapped)
					throw hash_table_overflow_exception();
				p = &table[0];
				wrapped = true;
			}
		}
		return p;
	}

	entry	*table;
	size_t	 size_;

};

template<typename _V> class PHash_table
{

public:

	struct entry
	{
		uint16_t finger_print;
		_V value;
	}; // __attribute__((packed));

	PHash_table() :
		size_(0)
	{}

	PHash_table(size_t size) :
		table(new entry[size]),
		size_(size)
	{
		memset(table.get(), 0, size_ * sizeof(entry));
	}

	PHash_table(size_t size, double factor) :
		size_(std::max(size_t((double)size * factor), (size_t)1llu)),
		table(new entry[size_])
	{
		memset(table.get(), 0, size_ * sizeof(entry));
	}
	
	entry* operator[](uint64_t hash) const
	{
		entry *entry = get_entry(hash);
		if (entry->value == 0u)
			return NULL;
		return entry;
	}

	entry* insert(uint64_t hash)
	{
		entry *entry = get_entry(hash);
		if (entry->value == 0u)
			entry->finger_print = finger_print(hash);
		return entry;
	}

	size_t size() const
	{
		return size_;
	}

	size_t count() const
	{
		size_t n(0);
		for (size_t i = 0; i<size_; ++i)
			if (table[i].value != 0u)
				++n;
		return n;
	}

private:

	static uint16_t finger_print(uint64_t hash)
	{
		return (uint16_t)(hash & 0xffffllu);
	}

	entry* get_entry(uint64_t hash) const
	{
		entry *p = table.get() + ((hash >> 16) % size_);
		const uint16_t fp = finger_print(hash);
		bool wrapped = false;
		while (p->finger_print != fp && p->value != 0u) {
			++p;
			if (p == table.get() + size_) {
				if (wrapped)
					throw hash_table_overflow_exception();
				p = table.get();
				wrapped = true;
			}
		}
		return p;
	}

	size_t size_;
	auto_ptr<entry> table;	

};


struct PHash_set
{

public:

	typedef uint8_t fp;

	PHash_set() :
		size_(0)
	{}

	PHash_set(size_t size) :
		table(new fp[size]),
		size_(size)
	{
		memset(table.get(), 0, size_ * sizeof(fp));
	}

	bool contains(uint64_t key) const
	{
		fp *entry = get_entry(key);
		return *entry != 0;
	}

	void insert(uint64_t key)
	{
		fp *entry = get_entry(key);
		if (*entry == (fp)0)
			*entry = finger_print(murmur_hash()(key));
	}

	size_t size() const
	{
		return size_;
	}

private:

	static fp finger_print(uint64_t hash)
	{
		return std::max((fp)(hash & ((1llu<<(sizeof(fp)*8))-1llu)), (fp)1);
	}

	fp* get_entry(uint64_t key) const
	{
		const uint64_t hash = murmur_hash()(key), f = finger_print(hash);
		fp *p = table.get() + ((hash >> sizeof(fp)*8) % size_);
		bool wrapped = false;
		while (*p != f && *p != (fp)0) {
			++p;
			if (p == table.get() + size_) {
				if (wrapped)
					throw hash_table_overflow_exception();
				p = table.get();
				wrapped = true;
			}
		}
		return p;
	}

	auto_ptr<fp> table;
	size_t size_;

};

#endif /* HASH_TABLE_H_ */
