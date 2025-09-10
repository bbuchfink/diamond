/****
DIAMOND protein aligner
Copyright (C) 2019-2024 Max Planck Society for the Advancement of Science e.V.

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

#include <memory>
#include <numeric>
#include "../dp.h"
#include "config.h"
#include "anchored.h"
#include "../score_profile.h"
#include "util/simd/dispatch.h"
#include "stats/score_matrix.h"
#include "align/def.h"

using std::list;
using std::vector;
using std::unique_ptr;
using std::pair;
using std::sort;
using std::runtime_error;
using std::accumulate;
using std::make_unique;

namespace DP { namespace BandedSwipe { namespace DISPATCH_ARCH {

struct TargetVector {
	//vector<DP::AnchoredSwipe::Target<int8_t>> int8;
	vector<DP::AnchoredSwipe::Target<int16_t>> int16;
};

struct Profiles {
	struct Reverse {};
	Profiles(Sequence seq, const int8_t* cbs, int64_t padding, const ScoreMatrix& matrix)
	{
		int16 = make_profile16(seq, cbs, padding, &matrix);
	}
	Profiles(const Profiles& p, Reverse):
		int16(p.int16.reverse())
	{}
	//LongScoreProfile<int8_t> int8;
	LongScoreProfile<int16_t> int16;
};

static Loc get_band(Loc qlen, Extension::Mode extension_mode) {
	//return (Extension::band(qlen, extension_mode) + 15) / 16 * 16;
	return config.sensitivity >= Sensitivity::ULTRA_SENSITIVE ? 160 : (config.sensitivity >= Sensitivity::MORE_SENSITIVE ? 96 : 32);
}

static void align_right(Sequence target_seq, bool reverse, Loc i, Loc j, Loc d_begin, Loc d_end, Score prefix_score, TargetVector& targets,
	int64_t target_idx,
	const LongScoreProfile<int16_t>* profile, const LongScoreProfile<int16_t>* profile_rev,
	const DP::AnchoredSwipe::Config& cfg)
{
	Loc qlen = cfg.query.length() - i, tlen = target_seq.length();
	//const int band_cap = std::max(std::min(qlen / 2, tlen / 2), 1);
	//const int band = std::min(std::max(get_band(cfg.query.length()), Loc((d_end - d_begin) * 0.15)), band_cap);
	//const int band = std::max(get_band(cfg.query.length()), Loc((d_end - d_begin) * 0.15));
	const int band = get_band(cfg.query.length(), cfg.extension_mode);
	d_begin -= band;
	d_end += band - 1;
	const Loc d0 = Geo::clip_diag(Geo::diag_sub_matrix(d_begin, i, j), qlen, tlen),
		d1 = Geo::clip_diag(Geo::diag_sub_matrix(d_end, i, j), qlen, tlen);	
	tlen = std::min(tlen, Geo::j(qlen - 1, d0) + 1);
	assert(tlen > 0);
	assert(d1 >= d0);
	
	const Sequence clipped_target = reverse ? target_seq.subseq(target_seq.length() - tlen, target_seq.length()) : target_seq.subseq(0, tlen);

	auto& t = targets.int16;
	t.emplace_back(clipped_target, d0, d1 + 1, i, qlen, target_idx, reverse);
	t.back().profile = profile;
	t.back().profile_rev = profile_rev;
}

static void align_left(Sequence target_seq, Loc i, Loc j, Loc d_begin, Loc d_end, Score suffix_score, TargetVector& targets,
	int64_t target_idx,
	const LongScoreProfile<int16_t>* profile, const LongScoreProfile<int16_t>* profile_rev,
	const DP::AnchoredSwipe::Config& cfg)
{
	const Loc qlen = cfg.query.length(), tlen = target_seq.length();
	const Loc ir = qlen - 1 - i, jr = tlen - 1 - j;
	align_right(target_seq.subseq(0, j + 1), true, ir, jr, Geo::rev_diag(d_end, qlen, tlen), Geo::rev_diag(d_begin, qlen, tlen), suffix_score, targets, target_idx, profile, profile_rev, cfg);
}

static void add_target(DpTarget& t, TargetVector& targets,
	int64_t& target_idx,
	const DP::AnchoredSwipe::Config& cfg)
{
	//t.anchor = make_clipped_anchor(t.anchor, cfg.query, cfg.query_cbs, t.seq);
	//if (t.anchor.score == 0)
		//return;
	if (t.extend_right(cfg.query.length())) {
		const Loc i = t.anchor.query_end(), j = t.anchor.subject_end();
		align_right(t.seq.subseq(j), false, i, j, t.anchor.d_min_right, t.anchor.d_max_right, t.anchor.prefix_score, targets, target_idx, t.prof, t.prof_reverse, cfg);
		++target_idx;
	}
	if (t.extend_left()) {
		const Score suffix_score = cfg.score_hint - t.anchor.prefix_score + t.anchor.score;
		align_left(t.seq, t.anchor.query_begin() - 1, t.anchor.subject_begin() - 1, t.anchor.d_min_left, t.anchor.d_max_left, suffix_score,
			targets, target_idx, t.prof, t.prof_reverse, cfg);
		++target_idx;
	}
}

static void swipe_threads(DP::AnchoredSwipe::Target<int16_t>* targets, int64_t count, const DP::AnchoredSwipe::Options& options, const DP::AnchoredSwipe::Config& cfg) {
	using Target = DP::AnchoredSwipe::Target<int16_t>;
	ThreadPool::TaskSet task_set(*cfg.thread_pool, 0);
	int64_t size = 0;
	Target* i0 = targets, *i1 = targets, *end = targets + count;
	while (i1 < end) {
		const auto n = std::min((ptrdiff_t)16, end - i1);
		size += accumulate(i1, i1 + n, (int64_t)0, [](int64_t n, const Target& t) {return n + t.gross_cells(); });
		i1 += n;
		if (size >= config.swipe_task_size) {
#if ARCH_ID == 2
			task_set.enqueue(DP::AnchoredSwipe::DISPATCH_ARCH::smith_waterman<::DISPATCH_ARCH::ScoreVector<int16_t, 0>>, i0, i1 - i0, options);
#endif
			cfg.stats.inc(Statistics::SWIPE_TASKS_TOTAL);
			cfg.stats.inc(Statistics::SWIPE_TASKS_ASYNC);
			i0 = i1;
			size = 0;
		}
	}
	if (task_set.total() == 0) {
		cfg.stats.inc(Statistics::SWIPE_TASKS_TOTAL);
#if ARCH_ID == 2
		DP::AnchoredSwipe::DISPATCH_ARCH::smith_waterman<::DISPATCH_ARCH::ScoreVector<int16_t, 0>>(i0, i1 - i0, options);
#endif
		return;
	}
	if (i1 - i0 > 0) {
		cfg.stats.inc(Statistics::SWIPE_TASKS_TOTAL);
		cfg.stats.inc(Statistics::SWIPE_TASKS_ASYNC);
#if ARCH_ID == 2
		task_set.enqueue(DP::AnchoredSwipe::DISPATCH_ARCH::smith_waterman<::DISPATCH_ARCH::ScoreVector<int16_t, 0>>, i0, i1 - i0, options);
#endif
	}
	task_set.run();
}

const ScoreMatrix& select_matrix(int qlen) {
	//if (qlen < 35) return pam30;
	//if (qlen <= 50)return pam70;
	//if (qlen <= 85)return blosum80;
	return score_matrix;
}

list<Hsp> anchored_swipe(Targets& targets, const DP::AnchoredSwipe::Config& cfg, std::pmr::monotonic_buffer_resource& pool) {
#if ARCH_ID != 2
	throw runtime_error("Anchored SWIPE requires at least AVX2 support");
#endif
	TaskTimer total;

	TargetVector target_vec;
	int64_t target_count = 0, target_len = 0;
	Loc max_target_len = 0;
	if(!cfg.target_profiles)
		for (int bin = 0; bin < DP::BINS; ++bin)
			for (const DpTarget& t : targets[bin]) {
				++target_count;
				target_len += t.seq.length();
				max_target_len = std::max(max_target_len, t.seq.length());
			}
	else
		for (int bin = 0; bin < DP::BINS; ++bin)
			target_count += targets[bin].size();

	TaskTimer timer;
	//target_vec.int8.reserve(target_count * 2);
	target_vec.int16.reserve(target_count * 2); // check for 16 bit overflow
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_ALLOC, timer.microseconds());

	timer.go();
	unique_ptr<Profiles> profiles, profiles_rev;
	vector<const int16_t*> prof_pointers, prof_pointers_rev;
	const ScoreMatrix& matrix = score_matrix; // select_matrix(cfg.query.length());

	if (!cfg.target_profiles) {
		profiles = make_unique<Profiles>(cfg.query, cfg.query_cbs, cfg.query.length() + max_target_len + 32, matrix);
		profiles_rev = make_unique<Profiles>(*profiles, Profiles::Reverse());
		prof_pointers = profiles->int16.pointers(0);
		prof_pointers_rev = profiles_rev->int16.pointers(0);
	}
	cfg.stats.inc(Statistics::TIME_PROFILE, timer.microseconds());	

	timer.go();
	int64_t target_idx = 0;
	for (int bin = 0; bin < DP::BINS; ++bin)
		for (DpTarget& t : targets[bin]) {
			add_target(t, target_vec, target_idx, cfg);
		}
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_ADD, timer.microseconds());

	auto& t = target_vec.int16;
	
	timer.go();
	//sort(target_vec.int8.begin(), target_vec.int8.end());
	sort(target_vec.int16.begin(), target_vec.int16.end());
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_SORT, timer.microseconds());

	DP::AnchoredSwipe::Stats stats;
	DP::AnchoredSwipe::Options options{ prof_pointers.empty() ? nullptr : prof_pointers.data(), prof_pointers_rev.empty() ? nullptr : prof_pointers_rev.data() };

	timer.go();
#ifdef __SSE4_1__
	//DP::AnchoredSwipe::DISPATCH_ARCH::smith_waterman<::DISPATCH_ARCH::ScoreVector<int8_t, 0>>(target_vec.int8.data(), target_vec.int8.size());
	//stats = DP::AnchoredSwipe::DISPATCH_ARCH::smith_waterman<::DISPATCH_ARCH::ScoreVector<int16_t, 0>>(target_vec.int16.data(), target_vec.int16.size(), options);
	swipe_threads(target_vec.int16.data(), target_vec.int16.size(), options, cfg);
	//cfg.stats.inc(Statistics::GROSS_DP_CELLS, stats.gross_cells);
	//cfg.stats.inc(Statistics::NET_DP_CELLS, stats.net_cells);
#endif
	cfg.stats.inc(Statistics::TIME_SW, timer.microseconds());

	timer.go();
	//sort(target_vec.int8.begin(), target_vec.int8.end(), DP::AnchoredSwipe::Target<int8_t>::cmp_target_idx);
	sort(target_vec.int16.begin(), target_vec.int16.end(), DP::AnchoredSwipe::Target<int16_t>::cmp_target_idx);
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_SORT, timer.microseconds());

	timer.go();
	auto target_it = target_vec.int16.cbegin();
	list<Hsp> out;
	DP::Targets recompute;
	std::pmr::list<Stats::TargetMatrix> matrices(&pool);
	Stats::Composition query_comp(Stats::composition(cfg.query));
	for (int bin = 0; bin < DP::BINS; ++bin)
		for (const DpTarget& t : targets[bin]) {
			if (t.anchor.score == 0)
				continue;
			int score = t.anchor.score, i0 = t.anchor.query_begin(), i1 = t.anchor.query_end(), j0 = t.anchor.subject_begin(), j1 = t.anchor.subject_end();
#ifdef DP_STAT
			int64_t gross_cells = 0, net_cells = 0;
#endif
			if (t.extend_right(cfg.query.length())) {
				score += target_it->score;
				i1 += target_it->query_end;
				j1 += target_it->target_end;
#ifdef DP_STAT
				auto c = target_it->cells();
				gross_cells += c.first;
				net_cells += c.second;
#endif
				++target_it;
			}
			if (t.extend_left()) {
				score += target_it->score;
				i0 -= target_it->query_end;
				j0 -= target_it->target_end;
#ifdef DP_STAT
				auto c = target_it->cells();
				gross_cells += c.first;
				net_cells += c.second;
#endif
				++target_it;
			}
			// filter for cover
			const double qcov = double(i1 - i0) / cfg.query.length() * 100.0, tcov = double(j1 - j0) / t.seq.length() * 100.0;
			if (config.query_or_target_cover > 0 || config.query_cover > 0 || config.subject_cover > 0) {				
				if (!cfg.recompute_adjusted && (std::max(qcov, tcov) < config.query_or_target_cover
					|| qcov < config.query_cover || tcov < config.subject_cover))
					continue;
			}

			/*if (std::max(i0, cfg.query.length() - i1) > config.uncov_cap || std::max(j0, t.seq.length() - j1) > config.uncov_cap)
				continue;*/
			
			const double evalue = matrix.evalue(score, cfg.query.length(), t.seq.length()); // consider adjusted matrix
			if(!cfg.recompute_adjusted && evalue > config.max_evalue)
				continue;
			
			if (cfg.recompute_adjusted) {
				matrices.emplace_back(query_comp, cfg.query.length(), Stats::CBS::MATRIX_ADJUST, t.seq, cfg.stats, pool, Stats::eUserSpecifiedRelEntropy);
				recompute[bin].push_back(DpTarget(t.seq, t.true_target_len, t.d_begin, t.d_end, t.target_idx,
					cfg.query.length(), &matrices.back(), DpTarget::CarryOver(), t.anchor));
			}
			else {
				out.emplace_back();
				out.back().score = score;
				out.back().evalue = evalue;
				out.back().bit_score = matrix.bitscore(score);
				out.back().swipe_target = t.target_idx;
				out.back().query_range = { i0, i1 };
				out.back().query_source_range = out.back().query_range;
				out.back().subject_range = { j0, j1 };
				out.back().subject_source_range = out.back().subject_range;
#ifdef DP_STAT
				out.back().reserved1 = (int)gross_cells; // stats.gross_cells;
				out.back().reserved2 = (int)net_cells; // stats.net_cells;
#endif
				out.back().approx_id = out.back().approx_id_percent(cfg.query, t.seq);
			}
		}
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_OUTPUT, timer.microseconds());
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE, total.microseconds());
	if (cfg.recompute_adjusted) {
		DP::Params params{ cfg.query, nullptr, Frame(0), cfg.query.length(), nullptr, DP::Flags::NONE, false,
		0, 0, HspValues::COORDS, cfg.stats, cfg.thread_pool };
		return DP::BandedSwipe::swipe(recompute, params);
	}
	return out;
}

}

DISPATCH_3(std::list<Hsp>, anchored_swipe, Targets&, targets, const DP::AnchoredSwipe::Config&, cfg, std::pmr::monotonic_buffer_resource&, pool)

}}