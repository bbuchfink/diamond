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

#include "seed_histogram.h"

seedp_range current_range;

Partitioned_histogram::Partitioned_histogram()
{ }

struct Histogram_callback
{
	Histogram_callback(size_t seqp, vector<shape_histogram> &data, const Seed_set *filter):
		filter(filter)
	{
		for (unsigned s = 0; s < shapes.count(); ++s)
			ptr.push_back(data[s][seqp].begin());
	}
	bool operator()(uint64_t seed, uint64_t pos, size_t shape)
	{
		if (filter == 0 || filter->contains(seed))
			++ptr[shape][seed_partition(seed)];
		return true;
	}
	void finish() const
	{
	}
	vector<unsigned*> ptr;
	const Seed_set *filter;
};

Partitioned_histogram::Partitioned_histogram(const Sequence_set &seqs, const Seed_set *filter) :
	data_(shapes.count()),
	p_(seqs.partition(config.threads_))
{
	for (unsigned s = 0; s < shapes.count(); ++s) {
		data_[s].resize(p_.size() - 1);
		memset(data_[s].data(), 0, (p_.size() - 1)*sizeof(unsigned)*Const::seedp);
	}
	Ptr_vector<Histogram_callback> cb;
	for (size_t i = 0; i < p_.size() - 1; ++i)
		cb.push_back(new Histogram_callback(i, data_, filter));
	seqs.enum_seeds(cb, p_, 0, shapes.count());
}

size_t Partitioned_histogram::max_chunk_size() const
{
	size_t max = 0;
	::partition<unsigned> p(Const::seedp, config.lowmem);
	for (unsigned shape = 0; shape < shapes.count(); ++shape)
		for (unsigned chunk = 0; chunk < p.parts; ++chunk)
			max = std::max(max, hst_size(data_[shape], seedp_range(p.getMin(chunk), p.getMax(chunk))));
	return max;
}

void Partitioned_histogram::build_seq_partition(const Sequence_set &seqs,
	const unsigned seqp,
	const size_t begin,
	const size_t end,
	vector<char> &buf,
	const Seed_set *filter)
{
	for (size_t i = begin; i < end; ++i) {

		assert(i < seqs.get_length());
		if (seqs[i].length() == 0)
			continue;
		Reduction::reduce_seq(seqs[i], buf);

		for (unsigned s = 0; s < shapes.count(); ++s) {
			Seed_iterator it(buf, shapes[s]);
			//Contiguous_seed_iterator<5, 4> it(buf);
			unsigned *ptr = data_[s][seqp].begin();
			uint64_t seed;
			while (it.good())
				if (it.get(seed, shapes[s]))
					if (filter == 0 || filter->contains(seed))
						++ptr[seed_partition(seed)];
		}

	}
}