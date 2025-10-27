/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <vector>

template<typename T>
struct PtrVector : public std::vector<T*>
{

	using vector = std::vector<T*>;

	PtrVector():
		vector()
	{}
	
	PtrVector(typename vector::size_type n):
		vector(n)
	{}

	T& operator[](typename vector::size_type n)
	{ return *vector::operator[](n); }

	const T& operator[](typename vector::size_type n) const
	{ return *vector::operator[](n); }

	T*& get(typename vector::size_type n)
	{
		return vector::operator[](n);
	}

	T& back()
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