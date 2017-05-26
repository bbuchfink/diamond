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

#include "seed_set.h"

struct Seed_set_callback
{
	Seed_set_callback(vector<bool> &data):
		data(&data)
	{}
	void operator()(uint64_t seed, uint64_t pos, uint64_t shape)
	{
		(*data)[seed] = true;
	}
	void finish()
	{}
	vector<bool> *data;
};

Seed_set::Seed_set(const Sequence_set &seqs):
	data_((size_t)pow(1llu<<Reduction::reduction.bit_size(), shapes[0].length_))
{
	if (!shapes[0].contiguous())
		throw std::runtime_error("Contiguous seed required.");
	vector<Seed_set_callback> v;
	v.push_back(Seed_set_callback(data_));
	seqs.enum_seeds(v, seqs.partition(1), 0, 1);
}