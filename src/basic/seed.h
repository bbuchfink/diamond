/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
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

#ifndef SEED_H_
#define SEED_H_

#include <stdint.h>
#include "const.h"
#include "../util/hash_function.h"
#include "config.h"
#include "value.h"
#include "score_matrix.h"

typedef uint64_t Packed_seed;

inline unsigned seed_partition(Packed_seed s)
{
	return (unsigned)(s & (Const::seedp-1));
}

inline unsigned seed_partition_offset(Packed_seed s)
{
	return (unsigned)(s >> Const::seedp_bits);
}

struct Hashed_seed
{
	Hashed_seed()
	{}
	explicit Hashed_seed(uint64_t seed):
		hash(murmur_hash()(seed))
	{}
	unsigned partition() const
	{
		return unsigned(hash&(p - 1));
	}
	uint64_t offset() const
	{
		return hash >> p_bits;
	}
	operator uint64_t() const
	{
		return hash;
	}
	enum {
		p_bits = 10, p = 1 << p_bits
	};
	uint64_t hash;
};

struct Seed
{
	Letter& operator[](unsigned i)
	{
		return data_[i];
	}
	friend std::ostream& operator<<(std::ostream &str, const Seed &s)
	{
		for (unsigned i = 0; i < config.seed_weight; ++i)
			str << value_traits.alphabet[(size_t)s.data_[i]];
		return str;
	}
	int score(const Seed &rhs) const
	{
		int s = 0;
		for (unsigned i = 0; i < config.seed_weight; ++i)
			s += score_matrix(data_[i], rhs.data_[i]);
		return s;
	}
	uint64_t packed() const
	{
		uint64_t s = 0;
		for (unsigned i = 0; i < config.seed_weight; ++i) {
			s *= 20;
			s += data_[i];			
		}
		return s;
	}
	void enum_neighborhood(int treshold, vector<Seed> &out);
private:
	void enum_neighborhood(unsigned pos, int treshold, vector<Seed>& out, int score);
	Letter data_[Const::max_seed_weight];
};

#endif /* SEED_H_ */
