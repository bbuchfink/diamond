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

#include "shape.h"
#include "sequence.h"

struct Seed_iterator
{
	Seed_iterator(vector<char> &seq, const shape &sh):
		ptr_ (seq.data()),
		end_ (ptr_ + seq.size() - sh.length_ + 1)
	{}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed, const shape &sh)
	{
		return sh.set_seed_reduced(seed, ptr_++);
	}
private:
	const char *ptr_, *end_;
};