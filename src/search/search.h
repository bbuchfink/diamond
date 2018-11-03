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

#ifndef SEARCH_H_
#define SEARCH_H_

#include <stddef.h>

void search_shape(unsigned sid, unsigned query_block);
bool use_single_indexed(double coverage, size_t query_letters, size_t ref_letters);

extern const double SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE;

#endif