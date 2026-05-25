/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <thread>
#include "data/block/block.h"
#include "search/seed_complexity.h"
#include "util/ptr_vector.h"
#include "search/seed_array/seed_iterator.h"
#include "basic/shape_config.h"
#include "flags.h"

template<bool SKIP_SEED_POSITIONS>
struct SkipSeedPositionImpl;

template<>
struct SkipSeedPositionImpl<true>
{
	static inline bool apply(const EnumCfg& cfg, int shape, const uint64_t pos)
	{
		return (*cfg.skip_seed_positions)[shape].atomic_load(pos);
	}
};

template<>
struct SkipSeedPositionImpl<false>
{
	static inline bool apply(const EnumCfg&, int, const uint64_t)
	{
		return false;
	}
};

template<bool SKIP_SEED_POSITIONS>
static inline bool skip_seed_position(const EnumCfg& cfg, int shape, const uint64_t pos)
{
	return SkipSeedPositionImpl<SKIP_SEED_POSITIONS>::apply(cfg, shape, pos);
}

template<typename F, typename Filter, bool SKIP_SEED_POSITIONS>
//Search::SeedStats FLATTEN enum_seeds(SequenceSet* seqs, F* f, unsigned begin, unsigned end, const Filter* filter, const EnumCfg& cfg)
Search::SeedStats enum_seeds(SequenceSet* seqs, F* f, unsigned begin, unsigned end, const Filter* filter, const EnumCfg& cfg)
{
	uint64_t key;
	Search::SeedStats stats;
	GrowableBuffer<Letter> buf(300);
	for (unsigned i = begin; i < end; ++i) {
		if (UNLIKELY(cfg.skip) && (*cfg.skip)[i / align_mode.query_contexts])
			continue;
		if(UNLIKELY(config.min_query_len > 0) && seqs->source_length(i) < config.min_query_len)
			continue;
		const Sequence seq = (*seqs)[i];
		const auto l = seq.length();
		Reduction::reduce_seq(seq, buf);
		for (int shape_id = cfg.shape_begin; shape_id < cfg.shape_end; ++shape_id) {
			const Shape& sh = shapes[shape_id];
			if (UNLIKELY(l < sh.length_)) continue;
			SeedIterator<const Letter*> it(buf.data(), buf.data() + l, sh);
			Letter* ptr = seqs->ptr(i);
			Loc j = 0;
			while (it.good()) {
				if (it.get(key, sh)) {
					const uint64_t pos = seqs->position(i, j);
					if (!skip_seed_position<SKIP_SEED_POSITIONS>(cfg, shape_id, pos) && filter->contains(key, shape_id))
						(*f)(key, pos, i, shape_id);
				}
				++j;
			}
		}
	}
	f->finish();
	return stats;
}

template<typename F, typename Filter, typename It, bool SKIP_SEED_POSITIONS>
Search::SeedStats enum_seeds_minimizer(SequenceSet* seqs, F* f, unsigned begin, unsigned end, const Filter* filter, const EnumCfg& cfg, Loc it_param)
{
	Search::SeedStats stats;
	GrowableBuffer<Letter> buf(300);
	for (unsigned i = begin; i < end; ++i) {
		if (cfg.skip && (*cfg.skip)[i / align_mode.query_contexts])
			continue;
		if (config.min_query_len > 0 && seqs->source_length(i) < config.min_query_len)
			continue;
		const Sequence seq = (*seqs)[i];
        //if (align_mode.mode != AlignMode::blastn)
		const auto l = seq.length();
		Reduction::reduce_seq(seq, buf);
        //else
            //buf = seq.copy();
		for (int shape_id = cfg.shape_begin; shape_id < cfg.shape_end; ++shape_id) {
			const Shape& sh = shapes[shape_id];
			if (l < sh.length_) continue;
			It it(buf.data(), buf.data() + l, sh, it_param);
			while (it.good()) {
				const uint64_t key = *it;
				const uint64_t pos = seqs->position(i, it.pos());
				if (!skip_seed_position<SKIP_SEED_POSITIONS>(cfg, shape_id, pos) && filter->contains(key, shape_id))
					(*f)(key, pos, i, shape_id);
				++it;
			}
		}
	}
	f->finish();
	return stats;
}

template<typename F, uint64_t BITS, typename Filter, bool SKIP_SEED_POSITIONS>
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
				const uint64_t pos = seqs->position(i, Loc(it.seq_ptr(sh) - seq.data()));
				if (!skip_seed_position<SKIP_SEED_POSITIONS>(cfg, shape_id, pos) && filter->contains(key, shape_id)) {
					if (!cfg.filter_low_complexity_seeds || Search::seed_is_complex(it.seq_ptr(sh), sh, cfg.seed_cut))
						(*f)(key, pos, i, shape_id);
					else if (cfg.mask_low_complexity_seeds)
						*it.seq_ptr(sh) |= SEED_MASK;
				}
				++it;
			}
		}
	}
	f->finish();
}

template<typename F, typename It, typename Filter, bool SKIP_SEED_POSITIONS>
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
			if (it.get(key)) {
				const uint64_t pos = seqs->position(i, j);
				if (!skip_seed_position<SKIP_SEED_POSITIONS>(cfg, 0, pos) && filter->contains(key, 0))
					if ((*f)(key, pos, i, 0) == false)
						return;
			}
			++j;
		}
	}
	f->finish();
}

template<typename F, typename Filter, typename IteratorFilter, bool SKIP_SEED_POSITIONS>
static void enum_seeds_worker(F* f, SequenceSet* seqs, const unsigned begin, const unsigned end, const Filter* filter, Search::SeedStats* stats, const EnumCfg* cfg)
{
	static const char* errmsg = "Unsupported contiguous seed.";
	if (cfg->code == SeedEncoding::CONTIGUOUS) {
		const uint64_t b = Reduction::get_reduction().bit_size(), l = shapes[cfg->shape_begin].length_;
		switch (l) {
		case 7:
			switch (b) {
			case 4:
				enum_seeds_contiguous<F, ContiguousSeedIterator<7, 4, IteratorFilter>, Filter, SKIP_SEED_POSITIONS>(seqs, f, begin, end, filter, *cfg);
				break;
			default:
				throw std::runtime_error(errmsg);
			}
			break;
		case 6:
			switch (b) {
			case 4:
				enum_seeds_contiguous<F, ContiguousSeedIterator<6, 4, IteratorFilter>, Filter, SKIP_SEED_POSITIONS>(seqs, f, begin, end, filter, *cfg);
				break;
			default:
				throw std::runtime_error(errmsg);
			}
			break;
		case 5:
			switch (b) {
			case 4:
				enum_seeds_contiguous<F, ContiguousSeedIterator<5, 4, IteratorFilter>, Filter, SKIP_SEED_POSITIONS>(seqs, f, begin, end, filter, *cfg);
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
			enum_seeds_hashed<F, 4, Filter, SKIP_SEED_POSITIONS>(seqs, f, begin, end, filter, *cfg);
			break;
		default:
			throw std::runtime_error("Unsupported reduction.");
		}
	}
	else if(cfg->minimizer_window > 0)
		*stats = enum_seeds_minimizer<F, Filter, MinimizerIterator<const Letter*>, SKIP_SEED_POSITIONS>(seqs, f, begin, end, filter, *cfg, cfg->minimizer_window);
	else if(cfg->sketch_size > 0)
		*stats = enum_seeds_minimizer<F, Filter, SketchIterator, SKIP_SEED_POSITIONS>(seqs, f, begin, end, filter, *cfg, cfg->sketch_size);
	else
		*stats = enum_seeds<F, Filter, SKIP_SEED_POSITIONS>(seqs, f, begin, end, filter, *cfg);
}

template <typename F, typename Filter, bool SKIP_SEED_POSITIONS>
static Search::SeedStats enum_seeds_impl(Block& seqs, PtrVector<F>& f, const Filter* filter, const EnumCfg& cfg)
{
	std::vector<std::thread> threads;
	threads.reserve(f.size());
	std::vector<Search::SeedStats> stats(f.size());
	for (unsigned i = 0; i < f.size(); ++i) {
		const unsigned begin = (unsigned)((*cfg.partition)[i]), end = (unsigned)((*cfg.partition)[i + 1]);
		if (cfg.filter_masked_seeds)
			threads.emplace_back(enum_seeds_worker<F, Filter, FilterMaskedSeeds, SKIP_SEED_POSITIONS>, &f[i], &seqs.seqs(), begin, end, filter, &stats[i], &cfg);
		else
			threads.emplace_back(enum_seeds_worker<F, Filter, void, SKIP_SEED_POSITIONS>, &f[i], &seqs.seqs(), begin, end, filter, &stats[i], &cfg);
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

template <typename F, typename Filter>
Search::SeedStats enum_seeds(Block& seqs, PtrVector<F>& f, const Filter* filter, const EnumCfg& cfg)
{
	if (cfg.soft_masking != MaskingAlgo::NONE)
		seqs.soft_mask(cfg.soft_masking);
	if (cfg.skip_seed_positions)
		return enum_seeds_impl<F, Filter, true>(seqs, f, filter, cfg);
	return enum_seeds_impl<F, Filter, false>(seqs, f, filter, cfg);
}