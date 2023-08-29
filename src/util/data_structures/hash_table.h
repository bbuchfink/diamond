/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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
#include <stdexcept>
#include <stdlib.h>

struct Modulo {
	size_t operator()(size_t offset, size_t size) const {
		return offset % size;
	}
};

struct NoModulo {
	size_t operator()(size_t offset, size_t size) const {
		return offset;
	}
};

template<typename _K, typename _V, typename _HashFunction, typename ModuloOp>
struct HashTable : private _HashFunction
{

	struct Entry
	{
		_K key;
		_V value;
		bool blank() const
		{
			return value == (_V)0;
		}
	};

	HashTable(size_t size, const _HashFunction &hash) :
		_HashFunction(hash),
		table((Entry*)calloc(size, sizeof(Entry))),
		size_(size),
		destroy_(true)
	{
	}

	HashTable(char* table, size_t size):
		table((Entry*)table),
		size_(size / sizeof(Entry)),
		destroy_(false)
	{
	}

	~HashTable()
	{
		if (destroy_)
			free(table);
	}

	_V& operator[](_K key)
	{
		Entry *e = get_entry(key);
		if (e->blank())
			e->key = key;
		return *e;
	}

	_V* find(_K key)
	{
		return get_present_entry(key);
	}

	Entry* find_entry(_K key)
	{
		return get_present_entry(key);
	}

	Entry* insert(_K key)
	{
		return get_or_insert_entry(key);
	}

	size_t size() const
	{
		return size_;
	}

	size_t count() const
	{
		size_t n = 0;
		for (size_t i = 0; i < size_; ++i)
			if (!table[i].blank())
				++n;
		return n;
	}

	Entry* data()
	{
		return table;
	}

private:

	Entry* get_entry(_K key, bool stat=false)
	{
		Entry *p = &table[ModuloOp()(_HashFunction::operator()(key), size_)];
		bool wrapped = false;
		//if(stat) ++probe_n;
		while (p->key != key && !p->blank()) {
			//if(stat) ++probe_l;
			++p;
			if (p == &table[size_]) {
				if (wrapped)
					throw std::runtime_error("Hash table overflow.");
				p = &table[0];
				wrapped = true;
			}
		}
		return p;
	}

	Entry* get_present_entry(_K key, bool stat = false)
	{
		Entry *p = &table[ModuloOp()(_HashFunction::operator()(key), size_)];
		bool wrapped = false;
		//if(stat) ++probe_n;
		while(true) {
			//if(stat) ++probe_l;
			if (p->blank())
				return nullptr;
			if (p->key == key)
				return p;
			++p;
			if (p == &table[size_]) {
				if (wrapped)
					throw std::runtime_error("Hash table overflow.");
				p = &table[0];
				wrapped = true;
			}
		}
		return p;
	}

	Entry* get_or_insert_entry(_K key, bool stat = false)
	{
		Entry *p = &table[ModuloOp()(_HashFunction::operator()(key), size_)];
		bool wrapped = false;
		//if(stat) ++probe_n;
		while (true) {
			//if(stat) ++probe_l;
			if (p->key == key)
				return p;
			if (p->blank()) {
				p->key = key;
				return p;
			}
			++p;
			if (p == &table[size_]) {
				if (wrapped)
					throw std::runtime_error("Hash table overflow.");
				p = &table[0];
				wrapped = true;
			}
		}
		return p;
	}

	Entry *table;
	size_t size_;
	bool destroy_;

};
