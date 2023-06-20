/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
struct Const
{

	enum {
		build_version = 162,
#ifdef SINGLE_THREADED
		seedp_bits = 0,
#else
		seedp_bits = 10,
#endif
		seedp = 1<<seedp_bits,
		max_seed_weight = 32,
		max_shapes = 64
	};

	static const char* version_string;
	static const char* program_name;

};

#define SIMPLE_SEARCH
static constexpr int MAX_CONTEXT = 6;
// #define FREQUENCY_MASKING
// #define ST_JOIN
// #define NO_COLLISION_FILTER
