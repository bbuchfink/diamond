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
struct MurmurHash {
	uint64_t operator()(uint64_t h) const
	{
		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdLL;
		h ^= h >> 33;
		h *= 0xc4ceb9fe1a85ec53LL;
		h ^= h >> 33;
		return h;
	}
};

