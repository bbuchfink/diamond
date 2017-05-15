/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef PTR_VECTOR_H_
#define PTR_VECTOR_H_

#include <vector>

using std::vector;

template<typename _t>
struct Ptr_vector : public vector<_t*>
{

	Ptr_vector():
		vector<_t*>()
	{}
	
	Ptr_vector(typename vector<_t*>::size_type n):
		vector<_t*>(n)
	{}

	_t& operator[](typename vector<_t*>::size_type n)
	{ return *vector<_t*>::operator[](n); }

	const _t& operator[](typename vector<_t*>::size_type n) const
	{ return *vector<_t*>::operator[](n); }

	_t*& get(typename vector<_t*>::size_type n)
	{
		return vector<_t*>::operator[](n);
	}

	typename vector<_t*>::iterator erase(typename vector<_t*>::iterator first, typename vector<_t*>::iterator last)
	{
		for (typename vector<_t*>::iterator i = first; i < last; ++i)
			delete *i;
		return vector<_t*>::erase(first, last);
	}

	~Ptr_vector()
	{
		for(typename vector<_t*>::iterator i=this->begin(); i!=this->end(); ++i)
			delete *i;
	}

};

#endif /* PTR_VECTOR_H_ */
