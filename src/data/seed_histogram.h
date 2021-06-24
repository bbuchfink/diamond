/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include <limits>
#include <array>
#include "../basic/seed.h"
#include "sequence_set.h"
#include "../basic/shape_config.h"
#include "../basic/seed_iterator.h"
#include "enum_seeds.h"

typedef std::vector<std::array<unsigned, Const::seedp>> shape_histogram;

struct SeedPartitionRange
{
	SeedPartitionRange():
		begin_ (0),
		end_ (0)
	{ }
	SeedPartitionRange(unsigned begin, unsigned end):
		begin_ (begin),
		end_ (end)
	{ }
	bool contains(unsigned i) const
	{ return i >= begin_ && i < end_; }
	unsigned begin() const
	{ return begin_; }
	unsigned end() const
	{ return end_; }
	bool lower(unsigned i) const
	{ return i < begin_; }
	bool lower_or_equal(unsigned i) const
	{ return i < end_; }
	unsigned size() const
	{
		return end_ - begin_;
	}
	static SeedPartitionRange all()
	{
		return SeedPartitionRange(0, Const::seedp);
	}
private:
	unsigned begin_, end_;
};

extern SeedPartitionRange current_range;

inline size_t partition_size(const shape_histogram &hst, size_t p)
{
	size_t s = 0;
	for (unsigned i = 0; i < hst.size(); ++i)
		s += hst[i][p];
	return s;
}

inline size_t hst_size(const shape_histogram &hst, const SeedPartitionRange &range)
{
	size_t s = 0;
	for(unsigned i=range.begin();i<range.end();++i)
		s += partition_size(hst, i);
	return s;
}

struct Partitioned_histogram
{

	Partitioned_histogram();
	
	template<typename _filter>
	Partitioned_histogram(SequenceSet &seqs, bool serial, const _filter *filter, bool hashed_seeds, const std::vector<bool>* skip) :
		data_(shapes.count()),
		p_(seqs.partition(config.threads_))
	{
		for (unsigned s = 0; s < shapes.count(); ++s) {
			data_[s].resize(p_.size() - 1);
			for (std::array<unsigned, Const::seedp>& h : data_[s])
				h.fill(0);
		}
		PtrVector<Callback> cb;
		for (size_t i = 0; i < p_.size() - 1; ++i)
			cb.push_back(new Callback(i, data_));
		if (serial)
			for (unsigned s = 0; s < shapes.count(); ++s)
				enum_seeds(&seqs, cb, p_, s, s + 1, filter, hashed_seeds, skip);
		else
			enum_seeds(&seqs, cb, p_, 0, shapes.count(), filter, hashed_seeds, skip);
	}

	const shape_histogram& get(unsigned sid) const
	{ return data_[sid]; }

	size_t max_chunk_size(size_t index_chunks) const;

	const vector<size_t>& partition() const
	{
		return p_;
	}

private:

	struct Callback
	{
		Callback(size_t seqp, vector<shape_histogram> &data)
		{
			for (unsigned s = 0; s < shapes.count(); ++s)
				ptr.push_back(data[s][seqp].data());
		}
		bool operator()(uint64_t seed, uint64_t pos, uint32_t block_id, size_t shape)
		{
			++ptr[shape][seed_partition(seed)];
			return true;
		}
		void finish() const
		{
		}
		vector<unsigned*> ptr;
	};

	vector<shape_histogram> data_;
	vector<size_t> p_;

};