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

#include "seed_histogram.h"
#include "../util/util.h"
#include "../util/algo/partition.h"

SeedPartitionRange current_range;

Partitioned_histogram::Partitioned_histogram()
{ }

size_t Partitioned_histogram::max_chunk_size(size_t index_chunks) const
{
	size_t max = 0;
	::Partition<unsigned> p(Const::seedp, index_chunks);
	for (unsigned shape = 0; shape < shapes.count(); ++shape)
		for (unsigned chunk = 0; chunk < p.parts; ++chunk)
			max = std::max(max, hst_size(data_[shape], SeedPartitionRange(p.begin(chunk), p.end(chunk))));
	return max;
}