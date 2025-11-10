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

#include "../util/simd/dispatch.h"
#include "../stage2.h"
#include "kernel.h"
#include "kernel_self.h"
#include "kernel_mutual_cov.h"
#include "kernel_lin.h"

namespace Search { namespace DISPATCH_ARCH {

typedef void Stage1KernelPackedLoc(const ::PackedLoc*, ::int32_t, const ::PackedLoc*, ::int32_t, ::Search::WorkSet&);
typedef void Stage1KernelPackedLocId(const ::PackedLocId*, ::int32_t, const ::PackedLocId*, ::int32_t, ::Search::WorkSet&);

static Stage1KernelPackedLocId* stage1_dispatch(const Search::Config* cfg, PackedLocId) {
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

static Stage1KernelPackedLoc* stage1_dispatch(const Search::Config* cfg, PackedLoc) {
	return config.lin_stage1 ? stage1_query_lin
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

}}

namespace Search {
	
DISPATCH_3V(run_stage1, JoinIterator<PackedLoc>&, it, Search::WorkSet*, work_set, const Search::Config*, cfg)
DISPATCH_3V(run_stage1, JoinIterator<PackedLocId>&, it, Search::WorkSet*, work_set, const Search::Config*, cfg)

}