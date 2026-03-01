/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#include "../util/simd/dispatch.h"
#include "../stage2.h"
#include "kernel.h"
#include "kernel_self.h"
#include "kernel_mutual_cov.h"
#include "kernel_lin.h"

namespace Search { namespace DISPATCH_ARCH {

typedef void Stage1KernelPackedLoc(const ::PackedLoc*, ::uint_fast32_t, const ::PackedLoc*, ::uint_fast32_t, ::Search::WorkSet&);
typedef void Stage1KernelPackedLocId(const ::PackedLocId*, ::uint_fast32_t, const ::PackedLocId*, ::uint_fast32_t, ::Search::WorkSet&);

static Stage1KernelPackedLocId* stage1_dispatch(const Search::Config* cfg, PackedLocId) {
	if (config.lin_stage1_combo)
		return stage1_longest_combo_lin;
	if (config.lin_stage1_query) {
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

static Stage1KernelPackedLoc* stage1_dispatch(const Search::Config* cfg, PackedLoc) {
	return config.lin_stage1_query ? stage1_query_lin
		: (cfg->lin_stage1_target ? stage1_target_lin<PackedLoc>
			: (config.self && cfg->current_ref_block == 0 ? stage1_self<PackedLoc> : stage1<PackedLoc>));
}

void run_stage1(JoinIterator<PackedLoc>& it, Search::WorkSet* work_set, const Search::Config* cfg) {
	auto kernel = stage1_dispatch(cfg, PackedLoc());
	for (; it; ++it) {
		work_set->stats.inc(Statistics::SEEDS_HIT);
		kernel(it.r->begin(), (int32_t)it.r->size(), it.s->begin(), (int32_t)it.s->size(), *work_set);
	}
}

void run_stage1(JoinIterator<PackedLocId>& it, Search::WorkSet* work_set, const Search::Config* cfg) {
	auto kernel = stage1_dispatch(cfg, PackedLocId());
	for (; it; ++it) {
		work_set->stats.inc(Statistics::SEEDS_HIT);
		kernel(it.r->begin(), (int32_t)it.r->size(), it.s->begin(), (int32_t)it.s->size(), *work_set);
	}
}

bool keep_target_id(const Search::Config& cfg) {
#ifdef HIT_KEEP_TARGET_ID
	return true;
#else
	return cfg.min_length_ratio != 0.0 || config.global_ranking_targets || (config.self && cfg.current_ref_block == 0) || config.lin_stage1_combo;
#endif
}

}}

namespace Search {
	
DISPATCH_3V(run_stage1, JoinIterator<PackedLoc>&, it, Search::WorkSet*, work_set, const Search::Config*, cfg)
DISPATCH_3V(run_stage1, JoinIterator<PackedLocId>&, it, Search::WorkSet*, work_set, const Search::Config*, cfg)
DISPATCH_1(bool, keep_target_id, const Search::Config&, cfg)

}