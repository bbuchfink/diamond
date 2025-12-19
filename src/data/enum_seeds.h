/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <thread>
#include "block/block.h"
#include "search/seed_complexity.h"
#include "util/ptr_vector.h"
#include "basic/seed_iterator.h"
#include "basic/shape_config.h"
#include "flags.h"

template<typename F, typename Filter>
Search::SeedStats enum_seeds(SequenceSet* seqs, F* f, unsigned begin, unsigned end, const Filter* filter, const EnumCfg& cfg)
{
	std::pmr::monotonic_buffer_resource pool;
	uint64_t key;
	Search::SeedStats stats;
	for (unsigned i = begin; i < end; ++i) {
		if (cfg.skip && (*cfg.skip)[i / align_mode.query_contexts])
			continue;
		if(config.min_query_len > 0 && seqs->source_length(i) < config.min_query_len)
			continue;
		const Sequence seq = (*seqs)[i];
		std::pmr::vector<Letter> buf = Reduction::reduce_seq(seq, pool);
		for (int shape_id = cfg.shape_begin; shape_id < cfg.shape_end; ++shape_id) {
			const Shape& sh = shapes[shape_id];
			if (seq.length() < sh.length_) continue;
			SeedIterator<typename std::pmr::vector<Letter>::const_iterator> it(buf.cbegin(), buf.cend(), sh);
			Letter* ptr = seqs->ptr(i);
			Loc j = 0;
			while (it.good()) {
				if (it.get(key, sh))
					if (filter->contains(key, shape_id))
						(*f)(key, seqs->position(i, j), i, shape_id);
				++j;
			}
		}
	}
	f->finish();
	return stats;
}

template<typename F, typename Filter, typename It>
Search::SeedStats enum_seeds_minimizer(SequenceSet* seqs, F* f, unsigned begin, unsigned end, const Filter* filter, const EnumCfg& cfg, Loc it_param)
{
	Search::SeedStats stats;
	std::pmr::monotonic_buffer_resource pool;
	for (unsigned i = begin; i < end; ++i) {
		if (cfg.skip && (*cfg.skip)[i / align_mode.query_contexts])
			continue;
		if (config.min_query_len > 0 && seqs->source_length(i) < config.min_query_len)
			continue;
		const Sequence seq = (*seqs)[i];
        //if (align_mode.mode != AlignMode::blastn)
        std::pmr::vector<Letter> buf = Reduction::reduce_seq(seq, pool);
        //else
            //buf = seq.copy();
		for (int shape_id = cfg.shape_begin; shape_id < cfg.shape_end; ++shape_id) {
			const Shape& sh = shapes[shape_id];
			if (seq.length() < sh.length_) continue;
			It it(buf.cbegin(), buf.cend(), sh, it_param);
			while (it.good()) {
				const uint64_t key = *it;
				if (filter->contains(key, shape_id))
					(*f)(key, seqs->position(i, it.pos()), i, shape_id);
				++it;
			}
		}
	}
	f->finish();
	return stats;
}

template<typename F, uint64_t BITS, typename Filter>
void enum_seeds_hashed(SequenceSet* seqs, F* f, unsigned begin, unsigned end, const Filter* filter, const EnumCfg& cfg)
{
	for (unsigned i = begin; i < end; ++i) {
		if (cfg.skip && (*cfg.skip)[i / align_mode.query_contexts])
			continue;
		if (config.min_query_len > 0 && seqs->source_length(i) < config.min_query_len)
			continue;
		const Sequence seq = (*seqs)[i];
		for (int shape_id = cfg.shape_begin; shape_id < cfg.shape_end; ++shape_id) {
			const Shape& sh = shapes[shape_id];
			if (seq.length() < sh.length_) continue;
			//const __m128i shape_mask = sh.long_mask_sse_;
			HashedSeedIterator<BITS> it(seqs->ptr(i), seqs->length(i), sh);
			while (it.good()) {
				const uint64_t key = *it;
				if (filter->contains(key, shape_id)) {
					if (!cfg.filter_low_complexity_seeds || Search::seed_is_complex(it.seq_ptr(sh), sh, cfg.seed_cut))
						(*f)(key, seqs->position(i, Loc(it.seq_ptr(sh) - seq.data())), i, shape_id);
					else if (cfg.mask_low_complexity_seeds)
						*it.seq_ptr(sh) |= SEED_MASK;
				}
				++it;
			}
		}
	}
	f->finish();
}

template<typename F, typename It, typename Filter>
void enum_seeds_contiguous(SequenceSet* seqs, F* f, unsigned begin, unsigned end, const Filter* filter, const EnumCfg& cfg)
{
	uint64_t key;
	for (unsigned i = begin; i < end; ++i) {
		if (cfg.skip && (*cfg.skip)[i / align_mode.query_contexts])
			continue;
		if (config.min_query_len > 0 && seqs->source_length(i) < config.min_query_len)
			continue;
		const Sequence seq = (*seqs)[i];
		if (seq.length() < It::length()) continue;
		It it(seq);
		Loc j = 0;
		while (it.good()) {
			if (it.get(key))
				if (filter->contains(key, 0))
					if ((*f)(key, seqs->position(i, j), i, 0) == false)
						return;
			++j;
		}
	}
	f->finish();
}

template<typename F, typename Filter, typename IteratorFilter>
static void enum_seeds_worker(F* f, SequenceSet* seqs, const unsigned begin, const unsigned end, const Filter* filter, Search::SeedStats* stats, const EnumCfg* cfg)
{
	using It = std::pmr::vector<Letter>::const_iterator;
	static const char* errmsg = "Unsupported contiguous seed.";
	if (cfg->code == SeedEncoding::CONTIGUOUS) {
		const uint64_t b = Reduction::get_reduction().bit_size(), l = shapes[cfg->shape_begin].length_;
		switch (l) {
		case 7:
			switch (b) {
			case 4:
				enum_seeds_contiguous<F, ContiguousSeedIterator<7, 4, IteratorFilter>, Filter>(seqs, f, begin, end, filter, *cfg);
				break;
			default:
				throw std::runtime_error(errmsg);
			}
			break;
		case 6:
			switch (b) {
			case 4:
				enum_seeds_contiguous<F, ContiguousSeedIterator<6, 4, IteratorFilter>, Filter>(seqs, f, begin, end, filter, *cfg);
				break;
			default:
				throw std::runtime_error(errmsg);
			}
			break;
		case 5:
			switch (b) {
			case 4:
				enum_seeds_contiguous<F, ContiguousSeedIterator<5, 4, IteratorFilter>, Filter>(seqs, f, begin, end, filter, *cfg);
				break;
			default:
				throw std::runtime_error(errmsg);
			}
			break;
		default:
			throw std::runtime_error(errmsg);
		}
	}
	else if (cfg->code == SeedEncoding::HASHED) {
		const uint64_t b = Reduction::get_reduction().bit_size();
		switch (b) {
		case 4:
			enum_seeds_hashed<F, 4, Filter>(seqs, f, begin, end, filter, *cfg);
			break;
		default:
			throw std::runtime_error("Unsupported reduction.");
		}
	}
	else if(cfg->minimizer_window > 0)
		*stats = enum_seeds_minimizer<F, Filter, MinimizerIterator<It>>(seqs, f, begin, end, filter, *cfg, cfg->minimizer_window);
	else if(cfg->sketch_size > 0)
		*stats = enum_seeds_minimizer<F, Filter, SketchIterator>(seqs, f, begin, end, filter, *cfg, cfg->sketch_size);
	else
		*stats = enum_seeds<F, Filter>(seqs, f, begin, end, filter, *cfg);
}

template <typename F, typename Filter>
Search::SeedStats enum_seeds(Block& seqs, PtrVector<F>& f, const Filter* filter, const EnumCfg& cfg)
{
	if (cfg.soft_masking != MaskingAlgo::NONE)
		seqs.soft_mask(cfg.soft_masking);
	std::vector<std::thread> threads;
	std::vector<Search::SeedStats> stats(f.size());
	for (unsigned i = 0; i < f.size(); ++i) {
		const unsigned begin = (unsigned)((*cfg.partition)[i]), end = (unsigned)((*cfg.partition)[i + 1]);
		if (cfg.filter_masked_seeds)
			threads.emplace_back(enum_seeds_worker<F, Filter, FilterMaskedSeeds>, &f[i], &seqs.seqs(), begin, end, filter, &stats[i], &cfg);
		else
			threads.emplace_back(enum_seeds_worker<F, Filter, void>, &f[i], &seqs.seqs(), begin, end, filter, &stats[i], &cfg);
	}
	for (auto& t : threads)
		t.join();
	for (size_t i = 1; i < f.size(); ++i) {
		stats[0].good_seed_positions += stats[i].good_seed_positions;
		stats[0].low_complexity_seeds += stats[i].low_complexity_seeds;
	}
	if (cfg.soft_masking != MaskingAlgo::NONE) {
		int l = 0;
		for (int i = cfg.shape_begin; i < cfg.shape_end; ++i)
			l = std::max(l, shapes[i].length_);
		seqs.remove_soft_masking((int)l, cfg.mask_seeds);
	}
	return stats[0];
}