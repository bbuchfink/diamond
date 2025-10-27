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
#include <list>
#include <vector>
#include "basic/sequence.h"
#include "basic/match.h"
#include "basic/statistics.h"
#include "basic/config.h"
#include "data/sequence_set.h"
#include "stats/cbs.h"
#include "flags.h"
#include "util/parallel/thread_pool.h"
#include "align/def.h"
#include "dp/score_profile.h"

struct DpTarget
{
	struct CarryOver {
		CarryOver() :
			i1(0), j1(0), ident(0), len(0)
		{}
		CarryOver(int i1, int j1, int ident, int len) :
			i1(i1), j1(j1), ident(ident), len(len)
		{}
		int i1, j1, ident, len;
	};
	enum { BLANK = -1, MIN_LETTERS = 3 };
	static Loc banded_cols(const Loc qlen, const Loc tlen, const Loc d_begin, const Loc d_end) {
		const Loc pos = std::max(d_end - 1, 0) - (d_end - 1);
		const Loc d0 = d_begin;
		const Loc j1 = std::min(qlen - 1 - d0, tlen - 1) + 1;
		return j1 - pos;
	}
	DpTarget():
		d_begin(),
		d_end(),
		cols(),
		target_idx(BLANK),
		matrix(nullptr)
	{
	}
	//DpTarget(const Sequence &seq, int true_target_len, int d_begin, int d_end, Interval chaining_target_range, Score chaining_score, BlockId target_idx, int qlen, const Stats::TargetMatrix* matrix = nullptr, const CarryOver& carry_over = CarryOver(), const Anchor& anchor = Anchor()) :
	DpTarget(const Sequence& seq, int true_target_len, int d_begin, int d_end, BlockId target_idx, int qlen, const Stats::TargetMatrix* matrix = nullptr, const CarryOver& carry_over = CarryOver(), const Anchor& anchor = Anchor()) :
		seq(seq),
		d_begin(d_begin),
		d_end(d_end),
		cols(banded_cols(qlen, seq.length(), d_begin, d_end)),
		true_target_len(true_target_len),
		//chaining_target_range(chaining_target_range),
		//chaining_score(chaining_score),
		target_idx(target_idx),
		carry_over(carry_over),
		matrix(matrix),
		anchor(anchor)
	{
	}
	DpTarget(const Sequence& seq, int true_target_len, BlockId target_idx, const Stats::TargetMatrix* matrix = nullptr, const CarryOver& carry_over = CarryOver()):
		seq(seq),
		d_begin(),
		d_end(),
		cols(),
		true_target_len(true_target_len),
		//chaining_score(0),
		target_idx(target_idx),
		carry_over(carry_over),
		matrix(matrix)
	{
	}
	DpTarget(const std::pair<const Letter*, int64_t> seq) :
		seq(seq.first, seq.second),
		d_begin(),
		d_end(),
		cols(),
		true_target_len((Loc)seq.second),
		target_idx(BLANK),
		matrix(nullptr)
	{
	}
	int left_i1() const
	{
		return std::max(d_end - 1, 0);
	}
	int band() const {
		return d_end - d_begin;
	}
	bool operator<(const DpTarget &x) const
	{
		const int i = left_i1(), j = x.left_i1(), b1 = band(), b2 = x.band(), bin_b1 = b1 / config.band_bin, bin_b2 = b2 / config.band_bin,
			t1 = cols, t2 = x.cols, bin_t1 = t1 / config.col_bin, bin_t2 = t2 / config.col_bin;
		return bin_b1 < bin_b2 || (bin_b1 == bin_b2 && (bin_t1 < bin_t2 || (bin_t1 == bin_t2 && i < j)));
		//return i < j || (i == j && (target_idx < x.target_idx || (target_idx == x.target_idx && d_begin < x.d_begin)));
	}
	bool blank() const {
		return target_idx == BLANK;
	}
	bool adjusted_matrix() const {
		return matrix != nullptr;
	}
	int matrix_scale() const {
		return adjusted_matrix() ? config.cbs_matrix_scale : 1;
	}
	int64_t cells(DP::Flags flags, Loc qlen) const {
		return flag_any(flags, DP::Flags::FULL_MATRIX) ? (int64_t)seq.length() * (int64_t)qlen
			: int64_t(d_end - d_begin) * (int64_t)cols;
	}
	bool extend_right(Loc qlen) const {
		//return anchor.query_end() < qlen && anchor.subject_end() < seq.length()
		return std::min(qlen - anchor.query_end(), seq.length() - anchor.subject_end()) >= MIN_LETTERS;
	}
	bool extend_left() const {
		//return anchor.query_begin() > 0 && anchor.subject_begin() > 0
		return std::min(anchor.query_begin(), anchor.subject_begin()) >= MIN_LETTERS;
	}
	friend std::ostream& operator<<(std::ostream& s, const DpTarget& t) {
		s << t.seq << '\t' << t.d_begin << '\t' << t.d_end << '\t' << t.cols << '\t' << t.true_target_len << '\t' << t.target_idx;
		return s;
	}
	Sequence seq;
	Loc d_begin, d_end, cols, true_target_len;
	//Interval chaining_target_range;
	//Score chaining_score;
	BlockId target_idx;
	CarryOver carry_over;
	const Stats::TargetMatrix* matrix;
	Anchor anchor;
	const LongScoreProfile<int16_t>* prof, *prof_reverse;
};

struct DpStat
{
	DpStat():
		gross_cells(0),
		net_cells(0)
	{}
	DpStat& operator+=(DpStat &x)
	{
		mtx_.lock();
		gross_cells += x.gross_cells;
		net_cells += x.net_cells;
		mtx_.unlock();
		return *this;
	}
	size_t gross_cells, net_cells;
private:
	std::mutex mtx_;
};

extern DpStat dp_stat;

namespace DP {

struct Params {
	const Sequence query;
	const char* query_id;
	const Frame frame;
	const int query_source_len;
	const int8_t* const composition_bias;
	const Flags flags;
	const bool reverse_targets;
	Loc target_max_len;
	int swipe_bin;
	HspValues v;
	Statistics& stat;
	ThreadPool* thread_pool;
};

enum { BINS = 6, SCORE_BINS = 3, ALGO_BINS = 2 };

struct Traceback {};
struct ScoreOnly {};

struct TargetVec {
	using const_iterator = std::vector<DpTarget>::const_iterator;
	using iterator = std::vector<DpTarget>::iterator;
	TargetVec() :
		max_len_(0)
	{}
	iterator begin() {
		return targets_.begin();
	}
	iterator end() {
		return targets_.end();
	}
	const_iterator begin() const {
		return targets_.begin();
	}
	const_iterator end() const {
		return targets_.end();
	}
	DpTarget& operator[](int64_t i) {
		return targets_[i];
	}
	const DpTarget& operator[](int64_t i) const {
		return targets_[i];
	}
	DpTarget& front() {
		return targets_.front();
	}
	const DpTarget& front() const {
		return targets_.front();
	}
	DpTarget& back() {
		return targets_.back();
	}
	int64_t size() const {
		return targets_.size();
	}
	void reserve(int64_t size) {
		targets_.reserve(size);
	}
	void push_back(const DpTarget& t) {
		targets_.push_back(t);
		max_len_ = std::max(max_len_, t.seq.length());
	}
	void push_back(const TargetVec& v) {
		targets_.insert(targets_.end(), v.targets_.begin(), v.targets_.end());
		max_len_ = std::max(max_len_, v.max_len_);
	}
	bool empty() const {
		return targets_.empty();
	}
	void clear() {
		targets_.clear();
		max_len_ = 0;
	}
	Loc max_len() const {
		return max_len_;
	}
	template<typename... T>
	void emplace_back(T&&... args) {
		targets_.emplace_back(std::forward<T>(args)...);
		max_len_ = std::max(max_len_, targets_.back().seq.length());
	}
private:
	std::vector<DpTarget> targets_;
	Loc max_len_;
};

using Targets = std::array<TargetVec, BINS>;

struct NoCBS {
	constexpr void* operator[](int i) const { return nullptr; }
};

namespace AnchoredSwipe {

struct Config {
	Sequence query;
	const int8_t* query_cbs;
	Score score_hint;
	Statistics& stats;
	ThreadPool* thread_pool;
	bool recompute_adjusted;
	Extension::Mode extension_mode;
	bool target_profiles;
};

}
	
namespace Swipe {

//DECL_DISPATCH(std::list<Hsp>, swipe, (const sequence &query, const sequence *subject_begin, const sequence *subject_end, int score_cutoff))

}

namespace BandedSwipe {

std::list<Hsp> swipe(const Targets& targets, Params& params);
std::list<Hsp> swipe_set(const SequenceSet::ConstIterator begin, const SequenceSet::ConstIterator end, Params& params);
int bin(HspValues v, int query_len, int score, int ungapped_score, const int64_t dp_size, unsigned score_width, const Loc mismatch_est);
std::list<Hsp> anchored_swipe(Targets& targets, const DP::AnchoredSwipe::Config& cfg, std::pmr::monotonic_buffer_resource& pool);

}

}

std::list<Hsp> banded_3frame_swipe(const TranslatedSequence& query, Strand strand, std::vector<DpTarget>::iterator target_begin, std::vector<DpTarget>::iterator target_end, DpStat& stat, bool score_only, bool parallel);