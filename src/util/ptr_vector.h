/****
Copyright (c) 2016, Benjamin Buchfink
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
