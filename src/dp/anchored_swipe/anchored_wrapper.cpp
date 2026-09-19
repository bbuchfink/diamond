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

#include <memory>
#include <numeric>
#include <cfloat>
#include <mutex>
#include "../dp.h"
#include "config.h"
#include "anchored.h"
#include "../long_profile/score_profile.h"
#include "util/simd/dispatch.h"
#include "stats/score_matrix.h"
#include "align/def.h"
#include "stats/stats.h"
#include "util/memory/mem_profile.h"

using std::list;
using std::vector;
using std::unique_ptr;
using std::pair;
using std::sort;
using std::accumulate;

namespace DP { namespace BandedSwipe { namespace DISPATCH_ARCH {

struct TargetVector {
	//vector<DP::AnchoredSwipe::Target<int8_t>> int8;
	vector<DP::AnchoredSwipe::Target<int16_t>> int16;
};

// Score profiles (MATRIX_ROW_SCORES off) or padded query letters (on) of an extension, for the directions that are extended.
struct Profiles {
	const LongScoreProfile<int16_t>* forward, *reverse;
	const Letter* query, *query_rev;
};

static Loc band_max(const DP::AnchoredSwipe::Config& cfg) {
	return cfg.max_diag_spread + 2 * cfg.band;
}

static void align_right(Loc query_len, Sequence target_seq, bool reverse, Loc i, Loc j, Loc d_begin, Loc d_end, Score prefix_score, TargetVector& targets,
	int64_t target_idx,
	const Profiles& profiles,
	const DP::AnchoredSwipe::Config& cfg)
{
	Loc qlen = query_len - i, tlen = target_seq.length();
	//const int band_cap = std::max(std::min(qlen / 2, tlen / 2), 1);
	//const int band = std::min(std::max(get_band(cfg.query.length()), Loc((d_end - d_begin) * 0.15)), band_cap);
	//const int band = std::max(get_band(cfg.query.length()), Loc((d_end - d_begin) * 0.15));
	const int band = cfg.band; // query_len?
	d_begin -= band;
	d_end += band - 1;
	const Loc d0 = Geo::clip_diag(Geo::diag_sub_matrix(d_begin, i, j), qlen, tlen),
		d1 = Geo::clip_diag(Geo::diag_sub_matrix(d_end, i, j), qlen, tlen);	
	tlen = std::min(tlen, Geo::j(qlen - 1, d0) + 1);
	assert(tlen > 0);
	assert(d1 >= d0);
	assert(d1 >= 0 && d1 + 1 - d0 <= band_max(cfg));
	
	const Sequence clipped_target = reverse ? target_seq.subseq(target_seq.length() - tlen, target_seq.length()) : target_seq.subseq(0, tlen);

	auto& t = targets.int16;
	t.emplace_back(clipped_target, d0, d1 + 1, i, qlen, target_idx, reverse);
	t.back().profile = profiles.forward;
	t.back().profile_rev = profiles.reverse;
	t.back().query = profiles.query;
	t.back().query_rev = profiles.query_rev;
}

static void align_left(Loc query_len, Sequence target_seq, Loc i, Loc j, Loc d_begin, Loc d_end, Score suffix_score, TargetVector& targets,
	int64_t target_idx,
	const Profiles& profiles,
	const DP::AnchoredSwipe::Config& cfg)
{
	const Loc qlen = query_len, tlen = target_seq.length();
	const Loc ir = qlen - 1 - i, jr = tlen - 1 - j;
	align_right(qlen, target_seq.subseq(0, j + 1), true, ir, jr, Geo::rev_diag(d_end, qlen, tlen), Geo::rev_diag(d_begin, qlen, tlen), suffix_score, targets, target_idx, profiles, cfg);
}

static bool extend_right(const ExtensionPipeline::Extension& e) {
	return std::min(e.query.length() - e.anchor.query_end(), e.target.length() - e.anchor.subject_end()) >= DpTarget::MIN_LETTERS;
}

static bool extend_left(const ExtensionPipeline::Extension& e) {
	return std::min(e.anchor.query_begin(), e.anchor.subject_begin()) >= DpTarget::MIN_LETTERS;
}

static void add_target(const ExtensionPipeline::Extension& e, const Profiles& profiles, TargetVector& targets,
	int64_t& target_idx,
	const DP::AnchoredSwipe::Config& cfg)
{
	const Loc qlen = e.query.length();
	const Anchor& anchor = e.anchor;
	if (extend_right(e)) {
		const Loc i = anchor.query_end(), j = anchor.subject_end();
		align_right(qlen, e.target.subseq(j), false, i, j, anchor.d_min_right, anchor.d_max_right, anchor.prefix_score, targets, target_idx, profiles, cfg);
		++target_idx;
	}
	if (extend_left(e)) {
		const Score suffix_score = cfg.score_hint - anchor.prefix_score + anchor.score;
		align_left(qlen, e.target, anchor.query_begin() - 1, anchor.subject_begin() - 1, anchor.d_min_left, anchor.d_max_left, suffix_score,
			targets, target_idx, profiles, cfg);
		++target_idx;
	}
}

// Writes seq framed by `padding` padding letters on each side to out, reversed if reverse is set, and returns the pointer to its first letter.
static const Letter* pad_query(Sequence seq, bool reverse, int64_t padding, Letter* out) {
	std::fill(out, out + padding, DP::AnchoredSwipe::PADDING_LETTER);
	Letter* begin = out + padding;
	if (reverse)
		std::reverse_copy(seq.data(), seq.end(), begin);
	else
		std::copy(seq.data(), seq.end(), begin);
	std::fill(begin + seq.length(), begin + seq.length() + padding, DP::AnchoredSwipe::PADDING_LETTER);
	return begin;
}

void anchored_swipe(ExtensionPipeline::ExtensionQueue& queue, const DP::AnchoredSwipe::Config& cfg, const ExtensionPipeline::ExtensionCallback& callback) {
	MEM_SCOPE("dp/anchored-swipe");
	TaskTimer total;

	vector<ExtensionPipeline::Extension> extensions;
	{
		std::lock_guard<std::mutex> lock(queue.mtx);
		extensions.swap(queue.extensions);
	}
	if (extensions.empty())
		return;

	TargetVector target_vec;

	TaskTimer timer;
	target_vec.int16.reserve(extensions.size() * 2); // check for 16 bit overflow
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_ALLOC, timer.microseconds());

	timer.go();
	// Each extension carries its own query, so the score profiles are built per extension, and only for
	// the directions that are extended. The padding covers every target within the band limits of the call.
	list<LongScoreProfile<int16_t>> profile_storage;
	vector<Letter> query_storage;
	vector<Profiles> profiles(extensions.size(), Profiles{ nullptr, nullptr, nullptr, nullptr });
	const ScoreMatrix& matrix = score_matrix;
	alignas(32) int8_t score_table[32 * 32];
	if (DP::AnchoredSwipe::MATRIX_ROW_SCORES) {
		MEM_SCOPE("dp/query-profile");
		std::copy(matrix.matrix8(), matrix.matrix8() + 32 * 32, score_table);
		for (int l = 0; l < 32; ++l)
			score_table[(l << 5) + DP::AnchoredSwipe::PADDING_LETTER] = DP::AnchoredSwipe::PADDING_SCORE;
		const int64_t padding = std::max(DP::AnchoredSwipe::DISPATCH_ARCH::profile_padding<::DISPATCH_ARCH::ScoreVector<int16_t, 0>>(band_max(cfg)),
			(int64_t)LongScoreProfile<int16_t>::DEFAULT_PADDING);
		size_t size = 0;
		for (const ExtensionPipeline::Extension& e : extensions)
			if (e.anchor.score > 0)
				size += (size_t(extend_right(e)) + size_t(extend_left(e))) * (e.query.length() + 2 * padding);
		query_storage.resize(size);
		Letter* out = query_storage.data();
		for (size_t k = 0; k < extensions.size(); ++k) {
			const ExtensionPipeline::Extension& e = extensions[k];
			if (e.anchor.score <= 0)
				continue;
			assert(e.diag_spread() <= cfg.max_diag_spread);
			if (extend_right(e)) {
				profiles[k].query = pad_query(e.query, false, padding, out);
				out += e.query.length() + 2 * padding;
			}
			if (extend_left(e)) {
				profiles[k].query_rev = pad_query(e.query, true, padding, out);
				out += e.query.length() + 2 * padding;
			}
		}
	}
	else {
		MEM_SCOPE("dp/query-profile");
		const int64_t padding = DP::AnchoredSwipe::DISPATCH_ARCH::profile_padding<::DISPATCH_ARCH::ScoreVector<int16_t, 0>>(band_max(cfg));
		for (size_t k = 0; k < extensions.size(); ++k) {
			const ExtensionPipeline::Extension& e = extensions[k];
			const bool right = extend_right(e), left = extend_left(e);
			if (e.anchor.score <= 0 || (!right && !left))
				continue;
			assert(e.diag_spread() <= cfg.max_diag_spread);
			profile_storage.push_back(make_profile16(e.query, nullptr, padding, &matrix));
			if (right)
				profiles[k].forward = &profile_storage.back();
			if (left) {
				if (right)
					profile_storage.push_back(profile_storage.back().reverse());
				else
					for (int l = 0; l < AMINO_ACID_COUNT; ++l)
						std::reverse(profile_storage.back().data[l].begin(), profile_storage.back().data[l].end());
				profiles[k].reverse = &profile_storage.back();
			}
		}
	}
	cfg.stats.inc(Statistics::TIME_PROFILE, timer.microseconds());

	timer.go();
	int64_t target_idx = 0;
	for (size_t k = 0; k < extensions.size(); ++k)
		if (extensions[k].anchor.score > 0)
			add_target(extensions[k], profiles[k], target_vec, target_idx, cfg);
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_ADD, timer.microseconds());

	timer.go();
	sort(target_vec.int16.begin(), target_vec.int16.end());
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_SORT, timer.microseconds());

	DP::AnchoredSwipe::Stats stats;
	DP::AnchoredSwipe::Options options{ nullptr, nullptr, score_table };

	timer.go();
	stats = DP::AnchoredSwipe::DISPATCH_ARCH::smith_waterman<::DISPATCH_ARCH::ScoreVector<int16_t, 0>>(target_vec.int16.data(), target_vec.int16.size(), band_max(cfg), options);
	cfg.stats.inc(Statistics::SWIPE_TASKS_TOTAL);
	cfg.stats.inc(Statistics::TIME_SW, timer.microseconds());

	timer.go();
	sort(target_vec.int16.begin(), target_vec.int16.end(), DP::AnchoredSwipe::Target<int16_t>::cmp_target_idx);
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_SORT, timer.microseconds());

	timer.go();
	auto target_it = target_vec.int16.cbegin();
	for (ExtensionPipeline::Extension& e : extensions) {
		e.score = 0;
		e.evalue = DBL_MAX;
		if (e.anchor.score <= 0)
			continue;
		const Sequence query = e.query;
		int score = e.anchor.score, i0 = e.anchor.query_begin(), i1 = e.anchor.query_end(), j0 = e.anchor.subject_begin(), j1 = e.anchor.subject_end();
		if (extend_right(e)) {
			score += target_it->score;
			i1 += target_it->query_end;
			j1 += target_it->target_end;
			++target_it;
		}
		if (extend_left(e)) {
			score += target_it->score;
			i0 -= target_it->query_end;
			j0 -= target_it->target_end;
			++target_it;
		}
		// filter for cover
		cfg.stats.inc(Statistics::EXT16);
		const double qcov = double(i1 - i0) / query.length() * 100.0, tcov = double(j1 - j0) / e.target.length() * 100.0;
		if (config.query_or_target_cover > 0 || config.query_cover > 0 || config.subject_cover > 0) {
			if (std::max(qcov, tcov) < config.query_or_target_cover || qcov < config.query_cover || tcov < config.subject_cover)
				continue;
		}

		const double evalue = matrix.evalue(score, query.length(), e.target.length()); // consider adjusted matrix
		if (evalue > config.max_evalue)
			continue;

		e.score = score;
		e.evalue = evalue;
		e.query_range = { i0, i1 };
		e.target_range = { j0, j1 };
	}
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_OUTPUT, timer.microseconds());
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE, total.microseconds());

	callback(extensions);
}

}

DISPATCH_3V(anchored_swipe, ExtensionPipeline::ExtensionQueue&, queue, const DP::AnchoredSwipe::Config&, cfg, const ExtensionPipeline::ExtensionCallback&, callback)

}}