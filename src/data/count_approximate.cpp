/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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