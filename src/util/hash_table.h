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

#ifndef HASH_TABLE_H_
#define HASH_TABLE_H_

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
	} __attribute__((packed));

	hash_table(size_t size):
		table (new entry[size]),
		size_ (size)
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
		if(_E()(entry->value))
			return NULL;
		return entry;
	}

	void insert(_K key, _V value)
	{
		entry *entry = get_entry(key);
		if(_E()(entry->value))
			entry->key = key;
		entry->value = value;
	}

	size_t size() const
	{
		return size_;
	}

	size_t count() const
	{
		size_t n (0);
		for(size_t i=0;i<size_;++i)
			if(!_E()(table[i].value))
				++n;
		return n;
	}

private:

	hash_table(const hash_table &):
		table (NULL),
		size_ (0)
	{
		assert(false);
	}

	hash_table& operator=(const hash_table &)
	{
		assert(false);
		return *this;
	}

	entry* get_entry(_K key) const
	{
		entry *p = &table[_H()(key) % size_];
		bool wrapped = false;
		while(p->key != key && !_E()(p->value)) {
			++p;
			if(p == &table[size_]) {
				if(wrapped)
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

#endif /* HASH_TABLE_H_ */
