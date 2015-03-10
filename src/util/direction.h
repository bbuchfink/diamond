/***
	DIAMOND protein sequence aligner
    Copyright (C) 2014 Benjamin Buchfink

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

#ifndef DIRECTION_H_
#define DIRECTION_H_

struct Right { };
struct Left { };

template<typename _val>
_val get_dir(const _val* x, int i, const Right&)
{ return *(x+i); }

template<typename _val>
_val get_dir(const _val* x, int i, const Left&)
{ return *(x-i); }

template<typename _val>
const _val* get_dir_ptr(const _val* x, int i, const Right&)
{ return x+i; }

template<typename _val>
const _val* get_dir_ptr(const _val* x, int i, const Left&)
{ return x-i; }

template<typename _val>
const _val* inc_dir(const _val* x, const Right&)
{ return x+1; }

template<typename _val>
const _val* inc_dir(const _val* x, const Left&)
{ return x-1; }

#endif /* DIRECTION_H_ */
