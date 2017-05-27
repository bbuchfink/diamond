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

#include "sorted_list.h"
#include "../util/ptr_vector.h"
#include "queries.h"

char* sorted_list::alloc_buffer(const Partitioned_histogram &hst)
{
	return new char[sizeof(entry) * hst.max_chunk_size()];
}

sorted_list::sorted_list()
{}

sorted_list::const_iterator sorted_list::get_partition_cbegin(unsigned p) const
{
	return const_iterator(cptr_begin(p), cptr_end(p));
}

sorted_list::iterator sorted_list::get_partition_begin(unsigned p) const
{
	return iterator(ptr_begin(p), ptr_end(p));
}

sorted_list::Random_access_iterator sorted_list::random_access(unsigned p, size_t offset) const
{
	return Random_access_iterator(cptr_begin(p) + offset, cptr_end(p));
}

sorted_list::entry* sorted_list::ptr_begin(unsigned i) const
{
	return &data_[limits_[i]];
}

sorted_list::entry* sorted_list::ptr_end(unsigned i) const
{
	return &data_[limits_[i + 1]];
}

const sorted_list::entry* sorted_list::cptr_begin(unsigned i) const
{
	return &data_[limits_[i]];
}

const sorted_list::entry* sorted_list::cptr_end(unsigned i) const
{
	return &data_[limits_[i + 1]];
}

sorted_list::Ptr_set sorted_list::build_iterators(const shape_histogram &hst) const
{
	Ptr_set iterators(hst.size());
	for (unsigned i = 0; i < Const::seedp; ++i)
		iterators[0][i] = ptr_begin(i);

	for (unsigned i = 1; i < hst.size(); ++i)
		for (unsigned j = 0; j < Const::seedp; ++j)
			iterators[i][j] = iterators[i - 1][j] + hst[i - 1][j];
	return iterators;
}