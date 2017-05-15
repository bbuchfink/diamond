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

#include <set>
#include "index.h"
#include "../basic/shape_config.h"

template<unsigned PartitionBits> struct Flajolet_Martin_counter
{
	static const unsigned N = 1 << PartitionBits;
	static const uint64_t MASK = N - 1;
	static const uint64_t HIGH = MASK << (64 - PartitionBits);
	static const double PHI;
	uint64_t buckets[N];

	Flajolet_Martin_counter()
	{
		clear();
	}

	void clear()
	{
		for (unsigned i = 0; i<N; ++i)
			buckets[i] = 0;
	}

	void add(uint64_t hash)
	{
		buckets[hash & MASK] |= 1LL << ctz(HIGH | (hash >> PartitionBits));
	}

	double get()
	{
		int n = 0;
		for (unsigned i = 0; i<N; ++i)
			n += ctz(~buckets[i]);
		return double(N) / PHI * pow(2, double(n) / N);
	}

};

template<unsigned PartitionBits> const double Flajolet_Martin_counter<PartitionBits>::PHI = 0.77351;

struct Exact_counter
{
	Exact_counter():
		seeds(shape_to-shape_from)
	{}
	void operator()(Hashed_seed seed, size_t pos, unsigned shape_id)
	{
		seeds[shape_id - shape_from][seed.partition()].insert(seed);
	}
	void finish()
	{}
	vector<Array<std::set<uint64_t>, Hashed_seed::p> > seeds;
};

struct Approximate_counter
{
	Approximate_counter() :
		data(shape_to - shape_from)
	{}
	void operator()(Hashed_seed seed, size_t pos, unsigned shape_id)
	{
		data[shape_id - shape_from].add(seed);
	}
	void finish()
	{}
	enum { counter_pbits = 8 };
	vector<Flajolet_Martin_counter<counter_pbits> > data;
};

vector<Array<unsigned, Hashed_seed::p> > count_exact(const Sequence_set &seqs)
{
	vector<Exact_counter> counters(config.threads_);
	seqs.enum_seeds(counters);
	vector<Array<unsigned, Hashed_seed::p> > out(shape_to - shape_from);
	memset(out.data(), 0, (shape_to - shape_from)*Hashed_seed::p*sizeof(unsigned));
	for (unsigned s = 0; s < shape_to - shape_from; ++s)
		for (unsigned p = 0; p < Hashed_seed::p; ++p)
			for (unsigned t = 0; t < config.threads_; ++t)
				out[s][p] += (unsigned)counters[t].seeds[s][p].size();
	return out;	
}

vector<size_t> count_approximate(const Sequence_set &seqs)
{
	vector<Approximate_counter> counters(config.threads_);
	seqs.enum_seeds(counters);
	vector<size_t> out(shape_to - shape_from);
	for (unsigned s = 0; s < shape_to - shape_from; ++s)
		for (unsigned t = 0; t < config.threads_; ++t)
			out[s] += (size_t)counters[t].data[s].get();
	return out;
}