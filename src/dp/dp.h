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

#pragma once
#include <list>
#include <vector>
#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../stats/hauser_correction.h"
#include "../basic/statistics.h"
#include "../basic/config.h"
#include "../data/sequence_set.h"
#include "../stats/cbs.h"
#include "flags.h"
#include "../util/parallel/thread_pool.h"

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
	{}
	DpTarget(const Sequence &seq, int true_target_len, int d_begin, int d_end, Interval chaining_target_range, Score chaining_score, BlockId target_idx, int qlen, const Stats::TargetMatrix* matrix = nullptr, const CarryOver& carry_over = CarryOver(), const Anchor& anchor = Anchor()) :
		seq(seq),
		d_begin(d_begin),
		d_end(d_end),
		cols(banded_cols(qlen, seq.length(), d_begin, d_end)),
		true_target_len(true_target_len),
		chaining_target_range(chaining_target_range),
		chaining_score(chaining_score),
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
		chaining_score(0),
		target_idx(target_idx),
		carry_over(carry_over),
		matrix(matrix)
	{}
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
	Interval chaining_target_range;
	Score chaining_score;
	BlockId target_idx;
	CarryOver carry_over;
	const Stats::TargetMatrix* matrix;
	Anchor anchor;
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
	HspValues v;
	Statistics& stat;
	ThreadPool* thread_pool;
};

enum { BINS = 6, SCORE_BINS = 3, ALGO_BINS = 2 };

struct Traceback {};
struct ScoreOnly {};

using Targets = std::array<std::vector<DpTarget>, BINS>;

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
};

}
	
namespace Swipe {

//DECL_DISPATCH(std::list<Hsp>, swipe, (const sequence &query, const sequence *subject_begin, const sequence *subject_end, int score_cutoff))

}

namespace BandedSwipe {

DECL_DISPATCH(std::list<Hsp>, swipe, (const Targets& targets, Params& params))
DECL_DISPATCH(std::list<Hsp>, swipe_set, (const SequenceSet::ConstIterator begin, const SequenceSet::ConstIterator end, Params& params))
DECL_DISPATCH(unsigned, bin, (HspValues v, int query_len, int score, int ungapped_score, const int64_t dp_size, unsigned score_width, const Loc mismatch_est))
DECL_DISPATCH(std::list<Hsp>, anchored_swipe, (Targets& targets, const DP::AnchoredSwipe::Config& cfg))

}

}

DECL_DISPATCH(std::list<Hsp>, banded_3frame_swipe, (const TranslatedSequence &query, Strand strand, std::vector<DpTarget>::iterator target_begin, std::vector<DpTarget>::iterator target_end, DpStat &stat, bool score_only, bool parallel))
