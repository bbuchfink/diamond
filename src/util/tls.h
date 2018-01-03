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

#ifndef TLS_H_
#define TLS_H_

#include <vector>

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
				ptr_ = new std::vector<Ptr_wrapper_base*>;
			ptr_->push_back(new Ptr_wrapper<_t>(ptr));
		}
		return *ptr;
	}
	static void clear()
	{
		if (ptr_ == 0)
			return;
		for (std::vector<Ptr_wrapper_base*>::iterator i = ptr_->begin(); i != ptr_->end(); ++i)
			delete *i;
		delete ptr_;
	}
private:
	static TLS_PTR std::vector<Ptr_wrapper_base*> *ptr_;
};

#endif