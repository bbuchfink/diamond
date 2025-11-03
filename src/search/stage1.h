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

#include "search.h"
#include "stage2.h"
#include "data/block/block.h"
#include "hamming/kernel.h"

using std::vector;

namespace Search { namespace DISPATCH_ARCH {
	
using Container = vector<FingerPrint, Util::Memory::AlignmentAllocator<FingerPrint, 16>>;

static void all_vs_all_self(const FingerPrint* a, uint32_t na, FlatArray<uint32_t>& out, unsigned hamming_filter_id) {
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e = a[i];
		out.next();
		for (uint32_t j = i + 1; j < na; ++j)
			if (e.match(a[j]) >= hamming_filter_id)
				out.push_back(j);
	}
}

static void all_vs_all_mutual_cov(const PackedLocId* q, const PackedLocId* s, const FingerPrint* a, uint32_t na, const FingerPrint* b, uint32_t nb, FlatArray<uint32_t>& out, unsigned hamming_filter_id, WorkSet& work_set) {
	uint32_t j0 = 0, j1 = 0;
	const double mlr = work_set.cfg.min_length_ratio;
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e = a[i];
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		out.next();
		while (j0 < nb) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j0].block_id);
			const double lr = (double)qlen / tlen;
			if (lr >= mlr)
				break;
			++j0;
		}
		j1 = std::max(j1, j0);
		for (uint32_t j = j0; j < j1; ++j) {
			if (e.match(b[j]) >= hamming_filter_id)
				out.push_back(j);
		}
		for (; j1 < nb; ++j1) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j1].block_id);
			const double lr = (double)tlen / qlen;
			if (lr < mlr)
				break;
			if (e.match(b[j1]) >= hamming_filter_id)
				out.push_back(j1);
		}
	}
}

static void all_vs_all_self_mutual_cov(const PackedLocId* q, const FingerPrint* a, uint32_t na, FlatArray<uint32_t>& out, unsigned hamming_filter_id, WorkSet& work_set) {
	const double mlr = work_set.cfg.min_length_ratio;
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e = a[i];
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		out.next();
		for (uint32_t j = i + 1; j < na; ++j) {
			const Loc tlen = work_set.cfg.query->seqs().length(q[j].block_id);
			if ((double)tlen / qlen < mlr)
				break;
			work_set.stats.inc(Statistics::SEED_HITS);
			if (e.match(a[j]) >= hamming_filter_id)
				out.push_back(j);
		}
	}
}

template<typename SeedLoc>
static void load_fps(const SeedLoc* p, size_t n, Container& v, const SequenceSet& seqs)
{
	v.clear();
	v.reserve(n);
	const SeedLoc* end = p + n;
	for (; p < end; ++p) {
		v.emplace_back(seqs.data(*p));
	}
}

template<typename SeedLoc>
static void FLATTEN stage1(const SeedLoc* q, int32_t nq, const SeedLoc* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, &vs = work_set.vs;
#endif
	
	const int32_t tile_size = config.tile_size;
	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, nq * ns);
	load_fps(q, nq, vq, work_set.cfg.query->seqs());
	const int32_t qs = (int32_t)vq.size(), ss = (int32_t)vs.size();
	for (int32_t i = 0; i < qs; i += tile_size) {
		for (int32_t j = 0; j < ss; j += tile_size) {
			work_set.hits.clear();
			all_vs_all(vq.data() + i, std::min(tile_size, qs - i), vs.data() + j, std::min(tile_size, ss - j), work_set.hits, work_set.cfg.hamming_filter_id);
			search_tile(work_set.hits, i, j, q, s, work_set);
		}
	}
}

static void FLATTEN stage1_query_lin(const PackedLoc* q, int32_t nq, const PackedLoc* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, &vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	load_fps(q, 1, vq, work_set.cfg.query->seqs());
	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, ns);
	
	const int32_t ss = (int32_t)vs.size();
	for (int32_t j = 0; j < ss; j += tile_size) {
		work_set.hits.clear();
		all_vs_all(vq.data(), 1, vs.data() + j, std::min(tile_size, ss - j), work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, 0, j, q, s, work_set);
	}
}

static void FLATTEN stage1_query_lin_ranked(const PackedLocId* q, int32_t nq, const PackedLocId* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, &vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	const int32_t ranking = work_set.kmer_ranking->highest_ranking(q, q + nq);
	load_fps(q + ranking, 1, vq, work_set.cfg.query->seqs());
	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, ns);

	const int32_t ss = (int32_t)vs.size();
	for (int32_t j = 0; j < ss; j += tile_size) {
		work_set.hits.clear();
		all_vs_all(vq.data(), 1, vs.data() + j, std::min(tile_size, ss - j), work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, ranking, j, q, s, work_set);
	}
}

template<typename SeedLoc>
static void FLATTEN stage1_target_lin(const SeedLoc* q, int32_t nq, const SeedLoc* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	load_fps(q, nq, vq, work_set.cfg.query->seqs());
	load_fps(s, 1, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, nq);

	const int32_t qs = (int32_t)vq.size();
	for (int32_t j = 0; j < qs; j += tile_size) {
		work_set.hits.clear();
		all_vs_all(vq.data() + j, std::min(tile_size, qs - j), vs.data(), 1, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, j, 0, q, s, work_set);
	}
}

template<typename SeedLoc>
static void FLATTEN stage1_self(const SeedLoc* q, int32_t nq, const SeedLoc* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vs;
#else
	Container& vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	load_fps(s, ns, vs, work_set.cfg.target->seqs());

	work_set.stats.inc(Statistics::SEED_HITS, ns*(ns - 1) / 2);
	const int32_t ss = (int32_t)vs.size();
	for (int32_t i = 0; i < ss; i += tile_size) {
		work_set.hits.clear();
		all_vs_all_self(vs.data() + i, std::min(tile_size, ss - i), work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, i, i, s, s, work_set);
		for (int32_t j = i + tile_size; j < ss; j += tile_size) {
			work_set.hits.clear();
			all_vs_all(vs.data() + i, std::min(tile_size, ss - i), vs.data() + j, std::min(tile_size, ss - j), work_set.hits, work_set.cfg.hamming_filter_id);
			search_tile(work_set.hits, i, j, s, s, work_set);
		}
	}
}

static void FLATTEN stage1_self_mutual_cov(const PackedLocId* q, int32_t nq, const PackedLocId* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vs;
#else
	Container& vs = work_set.vs;
#endif

	load_fps(s, ns, vs, work_set.cfg.target->seqs());

	//work_set.stats.inc(Statistics::SEED_HITS, ns * (ns - 1) / 2);
	const int32_t ss = (int32_t)vs.size();
	work_set.hits.clear();
	all_vs_all_self_mutual_cov(s, vs.data(), ss, work_set.hits, work_set.cfg.hamming_filter_id, work_set);
	search_tile(work_set.hits, 0, 0, s, s, work_set);
}

static void FLATTEN stage1_mutual_cov_query_lin(const PackedLocId* q, int32_t nq, const PackedLocId* s, int32_t ns, WorkSet& work_set)
{
	const double mlr = work_set.cfg.min_length_ratio;
	const bool self = config.self && work_set.cfg.current_ref_block == 0;
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	load_fps(q, nq, vq, work_set.cfg.query->seqs());

	const int32_t qs = (int32_t)vq.size(), ss = (int32_t)vs.size();
	int32_t j = 0;
	for (int32_t i = 0; i < qs;) {
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		int32_t j1 = j;
		for (; j1 < ss; ++j1) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j1].block_id);
			if ((double)tlen / qlen < mlr)
				break;
		}
		work_set.hits.clear();
		//all_vs_all(vq.data() + i, 1, vs.data() + j, j1 - j, work_set.hits, work_set.cfg.hamming_filter_id);
		const int32_t qpos = self ? i + (j1 - j) / 2 : i;
		all_vs_all(vq.data() + qpos, 1, vs.data() + j, j1 - j, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, qpos, j, q, s, work_set);
		j = j1;
		if (j == ss)
			break;
		const Loc tlen = work_set.cfg.target->seqs().length(s[j].block_id);
		for (; i < qs; ++i) {
			const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
			if ((double)tlen / qlen >= mlr)
				break;
		}
	}
}

static void FLATTEN stage1_mutual_cov_target_lin(const PackedLocId* q, int32_t nq, const PackedLocId* s, int32_t ns, WorkSet& work_set)
{
	const double mlr = work_set.cfg.min_length_ratio;
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	load_fps(q, nq, vq, work_set.cfg.query->seqs());

	const int32_t qs = (int32_t)vq.size(), ss = (int32_t)vs.size();
	int32_t i = 0;
	for (int32_t j = 0; j < ss;) {
		const Loc tlen = work_set.cfg.target->seqs().length(s[j].block_id);
		int32_t i1 = i;
		for (; i1 < qs; ++i1) {
			const Loc qlen = work_set.cfg.query->seqs().length(q[i1].block_id);
			if ((double)qlen / tlen < mlr)
				break;
		}
		work_set.hits.clear();
		all_vs_all(vq.data() + i, i1 - i, vs.data() + j, 1, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, i, j, q, s, work_set);
		i = i1;
		if (i == qs)
			break;
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		for (; j < ss; ++j) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j].block_id);
			if ((double)qlen / tlen >= mlr)
				break;
		}
	}
}

static void FLATTEN stage1_mutual_cov(const PackedLocId* q, int32_t nq, const PackedLocId* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, nq * ns);
	load_fps(q, nq, vq, work_set.cfg.query->seqs());
	const int32_t qs = (int32_t)vq.size(), ss = (int32_t)vs.size();
	for (int32_t i = 0; i < qs; i += tile_size) {
		for (int32_t j = 0; j < ss; j += tile_size) {
			work_set.hits.clear();
			all_vs_all_mutual_cov(q + i, s + j, vq.data() + i, std::min(tile_size, qs - i), vs.data() + j, std::min(tile_size, ss - j), work_set.hits, work_set.cfg.hamming_filter_id, work_set);
			search_tile(work_set.hits, i, j, q, s, work_set);
		}
	}
}

typedef void Stage1FnPackedLoc(const ::PackedLoc*, ::int32_t, const ::PackedLoc*, ::int32_t, ::Search::WorkSet&);
typedef void Stage1FnPackedLocId(const ::PackedLocId*, ::int32_t, const ::PackedLocId*, ::int32_t, ::Search::WorkSet&);

static Stage1FnPackedLocId* stage1_dispatch(const Search::Config* cfg, PackedLocId) {
	if (config.lin_stage1) {
		return cfg->min_length_ratio > 0.0 ? stage1_mutual_cov_query_lin : stage1_query_lin_ranked;
	}
	if (cfg->lin_stage1_target) {
		return cfg->min_length_ratio > 0.0 ? stage1_mutual_cov_target_lin : stage1_target_lin<PackedLocId>;
	}
	if (cfg->min_length_ratio > 0.0) {
		return config.self && cfg->current_ref_block == 0 ? stage1_self_mutual_cov : stage1_mutual_cov;
	}
	if (config.self && cfg->current_ref_block == 0) {
		return stage1_self;
	}
	return stage1;
}

static Stage1FnPackedLoc* stage1_dispatch(const Search::Config* cfg, PackedLoc) {
	return config.lin_stage1 ? stage1_query_lin
		: (cfg->lin_stage1_target ? stage1_target_lin<PackedLoc>
			: (config.self && cfg->current_ref_block == 0 ? stage1_self<PackedLoc> : stage1<PackedLoc>));
}

void run_stage1(JoinIterator<PackedLoc>& it, Search::WorkSet* work_set, const Search::Config* cfg) {
	auto f = stage1_dispatch(cfg, PackedLoc());
	for (; it; ++it)
		f(it.r->begin(), (int32_t)it.r->size(), it.s->begin(), (int32_t)it.s->size(), *work_set);
}

void run_stage1(JoinIterator<PackedLocId>& it, Search::WorkSet* work_set, const Search::Config* cfg) {
	auto f = stage1_dispatch(cfg, PackedLocId());
	for (; it; ++it)
		f(it.r->begin(), (int32_t)it.r->size(), it.s->begin(), (int32_t)it.s->size(), *work_set);
}

}}