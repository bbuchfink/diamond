/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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

#define NOMINMAX
#include "../lib/mio/mmap.hpp"
#include <array>
#include "seed_set.h"
#include "../util/ptr_vector.h"
#include "../util/math/integer.h"
#include "enum_seeds.h"
#include "../util/system/system.h"
#include "../util/io/output_file.h"
#include "../basic/shape_config.h"

static const size_t PADDING = 32;
static const double HASH_TABLE_FACTOR = 1.25;

using std::vector;
using std::endl;
using std::get;
using std::runtime_error;
using std::string;

NoFilter no_filter;

struct Seed_set_callback
{
	Seed_set_callback(vector<bool> &data, size_t max_coverage):
		coverage(0),
		max_coverage(max_coverage),
		data(&data)
	{}
	bool operator()(uint64_t seed, uint64_t pos, uint32_t block_id, uint64_t shape)
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

SeedSet::SeedSet(Block &seqs, double max_coverage, const std::vector<bool>* skip, const double seed_cut, const MaskingAlgo soft_masking):
	data_((size_t)pow(1llu<<Reduction::reduction.bit_size(), shapes[0].length_))
{
	if (!shapes[0].contiguous())
		throw std::runtime_error("Contiguous seed required.");
	PtrVector<Seed_set_callback> v;
	v.push_back(new Seed_set_callback(data_, size_t(max_coverage*pow(Reduction::reduction.size(), shapes[0].length_))));
	const auto p = seqs.seqs().partition(1);
	const EnumCfg cfg{ &p, 0, 1, SeedEncoding::CONTIGUOUS, skip, true, false, seed_cut, soft_masking, 0, false, false };
	enum_seeds(seqs, v, &no_filter, cfg);
	coverage_ = (double)v.back().coverage / pow(Reduction::reduction.size(), shapes[0].length_);
}

struct Hashed_seed_set_callback
{
	Hashed_seed_set_callback(PtrVector<HashSet<Modulo2, Identity>> &dst):
		dst(dst)
	{}
	bool operator()(uint64_t seed, uint64_t pos, uint32_t block_id, uint64_t shape)
	{
		dst[shape].insert(seed);
		return true;
	}
	void finish()
	{}
	PtrVector<HashSet<Modulo2, Identity> > &dst;
};

HashedSeedSet::HashedSeedSet(Block &seqs, const std::vector<bool>* skip, const double seed_cut, const MaskingAlgo soft_masking)
{
	for (size_t i = 0; i < shapes.count(); ++i)
		data_.push_back(new Table(next_power_of_2(seqs.seqs().letters() * HASH_TABLE_FACTOR)));
	PtrVector<Hashed_seed_set_callback> v;
	v.push_back(new Hashed_seed_set_callback(data_));
	const auto p = seqs.seqs().partition(1);
	const EnumCfg cfg{ &p, 0, shapes.count(), SeedEncoding::HASHED, skip, false, false, seed_cut, soft_masking, 0, false, false };
	enum_seeds(seqs, v, &no_filter, cfg);

	vector<size_t> sizes;
	for (size_t i = 0; i < shapes.count(); ++i)
		sizes.push_back(data_[i].load());
	data_.clear();

	for (size_t i = 0; i < shapes.count(); ++i)
		data_.push_back(new Table(next_power_of_2(sizes[i] * HASH_TABLE_FACTOR)));
	enum_seeds(seqs, v, &no_filter, cfg);

	for (size_t i = 0; i < shapes.count(); ++i) {
		data_[i].finish();
		log_stream << "Shape=" << i << " Hash_table_size=" << data_[i].size() << " load=" << (double)data_[i].load() / data_[i].size() << endl;
	}
}

HashedSeedSet::HashedSeedSet(const string& index_file):
	mmap_(new mio::mmap_source(index_file))
{
	if (mmap_->length() < SEED_INDEX_HEADER_SIZE)
		throw runtime_error("Invalid seed index file.");
	const char* buf = mmap_->data();
	if (*(uint64_t*)buf != SEED_INDEX_MAGIC_NUMBER)
		throw runtime_error("Invalid seed index file.");
	if (*(uint32_t*)(buf + 8) != SEED_INDEX_VERSION)
		throw runtime_error("Invalid seed index file version.");
	uint32_t shape_count = *(uint32_t*)(buf + 12);
	if (shape_count != shapes.count())
		throw runtime_error("Index has a different number of shapes.");

	const size_t* size_ptr = (const size_t*)(buf + SEED_INDEX_HEADER_SIZE);
	uint8_t* data_ptr = (uint8_t*)(buf + SEED_INDEX_HEADER_SIZE + sizeof(size_t) * shape_count);

	for (unsigned i = 0; i < shapes.count(); ++i) {
		data_.push_back(new Table(data_ptr, *size_ptr));
		log_stream << "MMAPED Shape=" << i << " Hash_table_size=" << data_[i].size() << " load=" << (double)data_[i].load() / data_[i].size() << endl;
		data_ptr += *size_ptr + Table::PADDING;
		++size_ptr;
	}
}

HashedSeedSet::~HashedSeedSet() {
}

size_t HashedSeedSet::max_table_size() const {
	return (*std::max_element(data_.begin(), data_.end(), [](Table* a, Table* b) { return a->size() < b->size(); }))->size();
}