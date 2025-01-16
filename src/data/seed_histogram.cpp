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

#include "seed_histogram.h"
#include "util/algo/partition.h"
#include "basic/shape_config.h"
#include "block/block.h"
#include "util/ptr_vector.h"
#include "seed_set.h"
#include "data/enum_seeds.h"

using std::vector;

SeedPartitionRange current_range;

SeedHistogram::SeedHistogram()
{ }

size_t SeedHistogram::max_chunk_size(const int index_chunks) const
{
	size_t max = 0;
	::Partition<int> p(seedp(), index_chunks);
	for (unsigned shape = 0; shape < shapes.count(); ++shape)
		for (int chunk = 0; chunk < p.parts; ++chunk)
			max = std::max(max, hst_size(data_[shape], SeedPartitionRange(p.begin(chunk), p.end(chunk))));
	return max;
}

template<typename Filter>
SeedHistogram::SeedHistogram(Block& seqs, bool serial, const Filter* filter, EnumCfg& enum_cfg, int seedp_bits) :
	p_(seqs.seqs().partition(config.threads_)),
	data_(shapes.count(), ShapeHistogram(p_.size() - 1, vector<unsigned>(seedp_count(seedp_bits), 0)))
{
	struct Callback {
		Callback(size_t seqp, int seedp_bits, vector<ShapeHistogram>& data):
			seedp_mask(::seedp_mask(seedp_bits))
		{
			for (unsigned s = 0; s < shapes.count(); ++s)
				ptr.push_back(data[s][seqp].data());
		}
		bool operator()(PackedSeed seed, uint64_t pos, uint32_t block_id, size_t shape)
		{
			++ptr[shape][seed_partition(seed, seedp_mask)];
			return true;
		}
		void finish() const
		{
		}
		const PackedSeed seedp_mask;
		vector<unsigned*> ptr;
	};

	PtrVector<Callback> cb;
	for (size_t i = 0; i < p_.size() - 1; ++i)
		cb.push_back(new Callback(i, seedp_bits, data_));
	enum_cfg.partition = &p_;
	if (serial)
		for (unsigned s = 0; s < shapes.count(); ++s) {
			enum_cfg.shape_begin = s;
			enum_cfg.shape_end = s + 1;
			enum_seeds(seqs, cb, filter, enum_cfg);
		}
	else {
		enum_cfg.shape_begin = 0;
		enum_cfg.shape_end = shapes.count();
		enum_seeds(seqs, cb, filter, enum_cfg);
	}
}

template SeedHistogram::SeedHistogram(Block&, bool, const NoFilter*, EnumCfg&, int);
template SeedHistogram::SeedHistogram(Block&, bool, const SeedSet*, EnumCfg&, int);
template SeedHistogram::SeedHistogram(Block&, bool, const HashedSeedSet*, EnumCfg&, int);