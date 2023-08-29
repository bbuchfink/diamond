/****
DIAMOND protein aligner
Copyright (C) 2019-2023 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#include "search.h"
#include "stage2.h"

using std::vector;

namespace Search { namespace DISPATCH_ARCH {
	
using Container = vector<FingerPrint, Util::Memory::AlignmentAllocator<FingerPrint, 16>>;

static void all_vs_all(const FingerPrint* a, uint32_t na, const FingerPrint* b, uint32_t nb, FlatArray<uint32_t>& out, unsigned hamming_filter_id) {
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e = a[i];
		out.next();
		for (uint32_t j = 0; j < nb; ++j)
			if (e.match(b[j]) >= hamming_filter_id)
				out.push_back(j);
	}
}

static void all_vs_all_self(const FingerPrint* a, uint32_t na, FlatArray<uint32_t>& out, unsigned hamming_filter_id) {
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e = a[i];
		out.next();
		for (uint32_t j = i + 1; j < na; ++j)
			if (e.match(a[j]) >= hamming_filter_id)
				out.push_back(j);
	}
}

static void load_fps(const SeedLoc* p, size_t n, Container& v, const SequenceSet& seqs)
{
	v.clear();
	v.reserve(n);
	const SeedLoc* end = p + n;
	for (; p < end; ++p) {
		v.emplace_back(seqs.data(*p));
	}
}

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

static void FLATTEN stage1_query_lin(const SeedLoc* q, int32_t nq, const SeedLoc* s, int32_t ns, WorkSet& work_set)
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

static void FLATTEN stage1_query_lin_ranked(const SeedLoc* q, int32_t nq, const SeedLoc* s, int32_t ns, WorkSet& work_set)
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

void run_stage1(JoinIterator<SeedLoc>& it, Search::WorkSet* work_set, const Search::Config* cfg) {
#ifdef KEEP_TARGET_ID
	auto f = config.lin_stage1 ? stage1_query_lin_ranked
		: (cfg->lin_stage1_target ? stage1_target_lin
			: (config.self && cfg->current_ref_block == 0 ? stage1_self : stage1));
#else
	auto f = config.lin_stage1 ? stage1_query_lin
		: (cfg->lin_stage1_target ? stage1_target_lin
			: (config.self && cfg->current_ref_block == 0 ? stage1_self : stage1));
#endif
	for (; it; ++it)
		f(it.r->begin(), (int32_t)it.r->size(), it.s->begin(), (int32_t)it.s->size(), *work_set);
}

}}