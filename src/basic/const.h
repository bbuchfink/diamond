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

#ifndef CONST_H_
#define CONST_H_

struct Const
{

	enum {
		build_version = 132,
		seedp_bits = 10,
		seedp = 1<<seedp_bits,
		max_seed_weight = 32,
		max_shapes = 16,
	};

	static const char* version_string;
	static const char* program_name;
	static const char* id_delimiters;

};

#define SIMPLE_SEARCH
static constexpr int MAX_CONTEXT = 6;
// #define FREQUENCY_MASKING
// #define ST_JOIN
// #define NO_COLLISION_FILTER

#endif /* CONST_H_ */
