/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

template<typename _t>
struct PtrVector : public std::vector<_t*>
{

	using vector = std::vector<_t*>;

	PtrVector():
		vector()
	{}
	
	PtrVector(typename vector::size_type n):
		vector(n)
	{}

	_t& operator[](typename vector::size_type n)
	{ return *vector::operator[](n); }

	const _t& operator[](typename vector::size_type n) const
	{ return *vector::operator[](n); }

	_t*& get(typename vector::size_type n)
	{
		return vector::operator[](n);
	}

	_t& back()
	{
		return *vector::back();
	}

	typename vector::iterator erase(typename vector::iterator first, typename vector::iterator last)
	{
		for (typename vector::iterator i = first; i < last; ++i)
			delete *i;
		return vector::erase(first, last);
	}

	void clear()
	{
		for (typename vector::iterator i = this->begin(); i != this->end(); ++i)
			delete *i;
		vector::clear();
	}

	~PtrVector()
	{
		clear();
	}

};

