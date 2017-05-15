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

#ifndef DIRECTION_H_
#define DIRECTION_H_

struct Right { enum { mult = 1 }; };
struct Left { enum { mult = -1 }; };

Letter get_dir(const Letter* x, int i, const Right&)
{ return *(x+i); }

Letter get_dir(const Letter* x, int i, const Left&)
{ return *(x-i); }

const Letter* get_dir_ptr(const Letter* x, int i, const Right&)
{ return x+i; }

const Letter* get_dir_ptr(const Letter* x, int i, const Left&)
{ return x-i; }

const Letter* inc_dir(const Letter* x, const Right&)
{ return x+1; }

const Letter* inc_dir(const Letter* x, const Left&)
{ return x-1; }

#endif /* DIRECTION_H_ */
