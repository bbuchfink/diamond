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

#include "seed_set.h"
#include "../util/ptr_vector.h"
#include "../util/math/integer.h"
#include "enum_seeds.h"

using std::endl;

No_filter no_filter;

struct Seed_set_callback
{
	Seed_set_callback(vector<bool> &data, size_t max_coverage):
		coverage(0),
		max_coverage(max_coverage),
		data(&data)
	{}
	bool operator()(uint64_t seed, uint64_t pos, uint64_t shape)
	{
		if ((*data)[seed] == false) {
			(*data)[seed] = true;
			++coverage;
			if (coverage > max_coverage)
				return false;
		}
		return true;
	}
	void finish()
	{}
	size_t coverage, max_coverage;
	vector<bool> *data;
};

Seed_set::Seed_set(const Sequence_set &seqs, double max_coverage):
	data_((size_t)pow(1llu<<Reduction::reduction.bit_size(), shapes[0].length_))
{
	if (!shapes[0].contiguous())
		throw std::runtime_error("Contiguous seed required.");
	PtrVector<Seed_set_callback> v;
	v.push_back(new Seed_set_callback(data_, size_t(max_coverage*pow(Reduction::reduction.size(), shapes[0].length_))));
	enum_seeds(&seqs, v, seqs.partition(1), 0, 1, &no_filter, true);
	coverage_ = (double)v.back().coverage / pow(Reduction::reduction.size(), shapes[0].length_);
}

struct Hashed_seed_set_callback
{
	Hashed_seed_set_callback(PtrVector<PHash_set<Modulo2, No_hash>> &dst):
		dst(dst)
	{}
	bool operator()(uint64_t seed, uint64_t pos, uint64_t shape)
	{
		dst[shape].insert(seed);
		return true;
	}
	void finish()
	{}
	PtrVector<PHash_set<Modulo2, No_hash> > &dst;
};

Hashed_seed_set::Hashed_seed_set(const Sequence_set &seqs)
{
	for (size_t i = 0; i < shapes.count(); ++i)
		data_.push_back(new PHash_set<Modulo2, No_hash>(next_power_of_2(seqs.letters()*1.25)));
		//data_.push_back(new PHash_set<void, Modulo2>(seqs.letters() * 1.25));
	PtrVector<Hashed_seed_set_callback> v;
	v.push_back(new Hashed_seed_set_callback(data_));
	enum_seeds(&seqs, v, seqs.partition(1), 0, shapes.count(), &no_filter);

	vector<size_t> sizes;
	for (size_t i = 0; i < shapes.count(); ++i)
		sizes.push_back(data_[i].load());
	data_.clear();

	for (size_t i = 0; i < shapes.count(); ++i)
		data_.push_back(new PHash_set<Modulo2, No_hash>(next_power_of_2(sizes[i] * 1.25)));
		//data_.push_back(new PHash_set<void, Modulo2>(seqs.letters() * 1.25));
	enum_seeds(&seqs, v, seqs.partition(1), 0, shapes.count(), &no_filter);

	for (size_t i = 0; i < shapes.count(); ++i)
		log_stream << "Shape=" << i << " Hash_table_size=" << data_[i].size() << " load=" << (double)data_[i].load()/data_[i].size() << endl;
}