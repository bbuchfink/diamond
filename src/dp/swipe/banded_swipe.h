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

#include <array>
#include <algorithm>
#include <utility>
#include <list>
#include <limits.h>
#include "../dp.h"
#include "swipe.h"
#include "target_iterator.h"
#include "../../util/math/integer.h"
#include "../../basic/config.h"
#include "../util/data_structures/range_partition.h"
#include "../../util/intrin.h"
#include "../../util/memory/alignment.h"

using std::list;
using std::pair;
using std::vector;
using std::array;

using namespace DISPATCH_ARCH;

namespace DP { namespace BandedSwipe {
namespace DISPATCH_ARCH {

template<typename _sv, typename _cbs>
Hsp traceback(_cbs bias_correction, const TracebackMatrix<_sv> &dp, const DpTarget &target, int d_begin, typename ScoreTraits<_sv>::Score max_score, double evalue, int max_col, int channel, int i0, int i1, int max_band_i, Void, Params& p)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	const int j0 = i1 - (target.d_end - 1), d1 = target.d_end;
	typename TracebackMatrix<_sv>::TracebackIterator it(dp.traceback(max_col + 1, i0 + max_col, j0 + max_col, (int)p.query.length(), channel, max_score));
	
	Hsp out(true);
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<_sv>::int_score(max_score);
	out.evalue = evalue;
	out.bit_score = score_matrix.bitscore(out.score);
	out.corrected_bit_score = score_matrix.bitscore_corrected(out.score, p.query.length(), target.true_target_len);
	out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));
	out.matrix = target.matrix;

	out.frame = p.frame.index();
	out.d_begin = target.d_begin;
	out.d_end = target.d_end;
	out.query_range.end_ = it.i + 1;
	out.subject_range.end_ = it.j + 1;
	
	while (it.score() > ScoreTraits<_sv>::zero_score()) {
		const Letter q = p.query[it.i], s = target.seq[it.j];
		const Score m = score_matrix(q, s), score = it.score();
		const Score m2 = add_cbs_scalar(m, bias_correction[it.i]);
		if (score == saturated_add(it.diag(), m2)) {
			out.push_match(q, s, m > (Score)0);
			it.walk_diagonal();
		} else {
			const pair<Edit_operation, int> g(it.walk_gap(d_begin, d1));
			out.push_gap(g.first, g.second, target.seq.data() + it.j + g.second);
		}
	}

	out.query_range.begin_ = it.i + 1;
	out.subject_range.begin_ = it.j + 1;
	out.transcript.reverse();
	out.transcript.push_terminator();
	out.target_seq = target.seq;
	out.query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(out.query_range.begin_, p.frame), TranslatedPosition(out.query_range.end_, p.frame), p.query_source_len);
	out.approx_id = out.approx_id_percent(p.query, target.seq);
	return out;
}

template<typename _sv, typename Cell, typename _cbs, typename StatType>
Hsp traceback(_cbs bias_correction, const Matrix<Cell> &dp, const DpTarget &target, int d_begin, typename ScoreTraits<_sv>::Score max_score, double evalue, int max_col, int channel, int i0, int i1, int max_band_i, const StatType& stats, Params& p)
{
	Hsp out(false);
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<_sv>::int_score(max_score);
	if (!target.adjusted_matrix())
		out.score *= config.cbs_matrix_scale;
	out.evalue = evalue;
	out.bit_score = score_matrix.bitscore(out.score);
	out.corrected_bit_score = score_matrix.bitscore_corrected(out.score, p.query.length(), target.true_target_len);
	out.frame = p.frame.index();
	out.matrix = target.matrix;
	const int j0 = i1 - (target.d_end - 1), i1_ = i0 + max_col + max_band_i + 1, j1_ = j0 + max_col + 1;
	if (target.carry_over.i1 == 0) {
		out.d_begin = target.d_begin;
		out.d_end = target.d_end;
		out.query_range.end_ = i1_;
		out.subject_range.end_ = j1_;
		out.target_seq = target.seq;
	}
	else {
		out.d_begin = -target.d_end + (int)p.query.length() - (int)target.seq.length() + 1;
		out.d_end = -target.d_begin + (int)p.query.length() - (int)target.seq.length() + 1;
		out.query_range.end_ = target.carry_over.i1;
		out.subject_range.end_ = target.carry_over.j1;
		out.identities = target.carry_over.ident;
		out.length = target.carry_over.len;
		out.query_range.begin_ = (int)p.query.length() - i1_;
		out.subject_range.begin_ = (int)target.seq.length() - j1_;
		out.approx_id = out.approx_id_percent(Sequence(p.query.reverse()), Sequence(target.seq.reverse()));
	}
	assign_stats(out, stats);
	out.query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(out.query_range.begin_, p.frame), TranslatedPosition(out.query_range.end_, p.frame), p.query_source_len);
	return out;
}

template<typename _sv, typename _cbs>
Hsp traceback(_cbs bias_correction, const TracebackVectorMatrix<_sv> &dp, const DpTarget &target, int d_begin, typename ScoreTraits<_sv>::Score max_score, double evalue, int max_col, int channel, int i0, int i1, int max_band_i, Void, Params& p)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	typedef typename ScoreTraits<_sv>::TraceMask TraceMask;
	const auto channel_mask = TraceMask::vmask(channel) | TraceMask::hmask(channel);
	const int j0 = i1 - (target.d_end - 1), d1 = target.d_end;
	typename TracebackVectorMatrix<_sv>::TracebackIterator it(dp.traceback(max_col + 1, i0 + max_col, max_band_i, j0 + max_col, (int)p.query.length(), channel));
	Hsp out(true);
	out.swipe_target = target.target_idx;
	out.target_seq = target.seq;
	out.score = ScoreTraits<_sv>::int_score(max_score);
	out.evalue = evalue;
	out.bit_score = score_matrix.bitscore(out.score);
	out.corrected_bit_score = score_matrix.bitscore_corrected(out.score, p.query.length(), target.true_target_len);
	out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));
	out.matrix = target.matrix;

	out.frame = p.frame.index();
	out.d_begin = target.d_begin;
	out.d_end = target.d_end;
	out.query_range.end_ = it.i + 1;
	out.subject_range.end_ = it.j + 1;
	const int end_score = out.score;
	const bool adjusted_matrix = target.adjusted_matrix();
	if (!adjusted_matrix)
		out.score *= config.cbs_matrix_scale;
	int score = 0;
	const int* matrix = adjusted_matrix ? target.matrix->scores32.data() : score_matrix.matrix32();

	while (it.i >= 0 && it.j >= 0 && score < end_score) {
		if ((it.mask().gap & channel_mask) == 0) {
			const Letter q = p.query[it.i], s = target.seq[it.j];
			const int m = matrix[int(s) * 32 + (int)q];
			const int m2 = adjusted_matrix ? m : add_cbs_scalar(m, bias_correction[it.i]);
			score += m2;
			out.push_match(q, s, m > (Score)0);
			it.walk_diagonal();
		}
		else {
			const pair<Edit_operation, int> g(it.walk_gap());
			out.push_gap(g.first, g.second, target.seq.data() + it.j + g.second);
			score -= (score_matrix.gap_open() + g.second * score_matrix.gap_extend()) * target.matrix_scale();
		}
	}

	if (score != end_score)
		throw std::runtime_error("Traceback error.");

	out.query_range.begin_ = it.i + 1;
	out.subject_range.begin_ = it.j + 1;
	out.transcript.reverse();
	out.transcript.push_terminator();
	out.query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(out.query_range.begin_, p.frame), TranslatedPosition(out.query_range.end_, p.frame), p.query_source_len);
	out.approx_id = out.approx_id_percent(p.query, target.seq);
	return out;
}
template<typename _sv, typename _cbs, typename Cfg>
list<Hsp> swipe(const vector<DpTarget>::const_iterator subject_begin, const vector<DpTarget>::const_iterator subject_end, _cbs composition_bias, vector<DpTarget> &overflow, Params& p)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	using Cell = typename Cfg::Cell;
	using Matrix = typename SelectMatrix<Cell, Cfg::traceback>::Type;
	using RowCounter = typename Cfg::RowCounter;
	using IdMask = typename Cfg::IdMask;
	using StatType = decltype(extract_stats<_sv>(Cell(), int()));
 	constexpr int CHANNELS = ScoreTraits<_sv>::CHANNELS;

	assert(subject_end - subject_begin <= CHANNELS);
	const int qlen = (int)p.query.length();

	int band = 0;
	for (vector<DpTarget>::const_iterator j = subject_begin; j < subject_end; ++j)
		band = std::max(band, j->d_end - j->d_begin);

	if (band > RowCounter::MAX_LEN)
		throw std::runtime_error("Band size exceeds row counter maximum.");

	int i1 = INT_MAX, d_begin[CHANNELS];
	const int target_count = int(subject_end - subject_begin);
#ifdef STRICT_BAND
	int band_offset[CHANNELS];
#endif
	for (int i = 0; i < target_count; ++i) {
		d_begin[i] = subject_begin[i].d_end - band;
#ifdef STRICT_BAND
		band_offset[i] = subject_begin[i].d_begin - d_begin[i];
#endif		
		i1 = std::min(i1, std::max(subject_begin[i].d_end - 1, 0));
	}
	int i0 = i1 + 1 - band;
#ifdef STRICT_BAND
	RangePartition<CHANNELS, Score> band_parts(band_offset, target_count, band);
#endif
	
	::DISPATCH_ARCH::TargetIterator<Score> targets(subject_begin, subject_end, i1, qlen, d_begin);
	Matrix dp(band, targets.cols);

	const uint32_t cbs_mask = targets.cbs_mask();
	const Score go = score_matrix.gap_open() + score_matrix.gap_extend(), go_s = go * (Score)config.cbs_matrix_scale,
		ge = score_matrix.gap_extend(), ge_s = ge * (Score)config.cbs_matrix_scale;
	const _sv open_penalty = blend_sv<_sv>(go, go_s, cbs_mask),
		extend_penalty = blend_sv<_sv>(ge, ge_s, cbs_mask);
	SwipeProfile<_sv> profile;
	array<const int8_t*, 32> target_scores;

	Score best[CHANNELS];
	int max_col[CHANNELS], max_band_row[CHANNELS];
	StatType stats[CHANNELS];
	std::fill(best, best + CHANNELS, ScoreTraits<_sv>::zero_score());
	std::fill(max_col, max_col + CHANNELS, 0);
	std::fill(max_band_row, max_band_row + CHANNELS, 0);
	CBSBuffer<_sv, _cbs> cbs_buf(composition_bias, qlen, cbs_mask);

	int j = 0;
	while (targets.active.size() > 0) {
		const int i0_ = std::max(i0, 0), i1_ = std::min(i1, qlen - 1) + 1, band_offset = i0_ - i0;
		if (i0_ >= i1_)
			break;
		typename Matrix::ColumnIterator it(dp.begin(band_offset, j));
		Cell vgap = Cell(), hgap = Cell();
		_sv col_best = _sv();
		RowCounter row_counter(band_offset);

		if (band_offset > 0)
			it.set_zero();

		const auto target_seqv = targets.get();
		const _sv target_seq = _sv(target_seqv);
		if (cbs_mask != 0) {
			if (targets.custom_matrix_16bit)
				profile.set(targets.get32().data());
			else
				profile.set(targets.get(target_scores.data()));
		}
		else {
#ifdef __SSSE3__
			profile.set(target_seqv);
#else
			profile.set(targets.get(target_scores.data()));
#endif
		}
#ifdef DP_STAT
		const uint64_t live = targets.live();
#endif

#ifdef STRICT_BAND
		for (int part = 0; part < band_parts.count(); ++part) {
			const int i_begin = std::max(i0 + band_parts.begin(part), i0_);
			const int i_end = std::min(i0 + band_parts.end(part), i1_);
			const _sv target_mask = load_sv<_sv>(band_parts.mask(part));
			vgap += target_mask;
#ifdef DP_STAT
			p.stat.inc(Statistics::GROSS_DP_CELLS, uint64_t(i_end - i_begin) * CHANNELS);
			p.stat.inc(Statistics::NET_DP_CELLS, uint64_t(i_end - i_begin) * popcount64(live & band_parts.bit_mask(part)));
#endif
			for (int i = i_begin; i < i_end; ++i) {
#else
			for (int i = i0_; i < i1_; ++i) {
#endif
				hgap = it.hgap();
				auto stat_h = it.hstat();
				_sv match_scores = profile.get(p.query[i]);
#ifdef STRICT_BAND
				hgap += target_mask;
				match_scores += target_mask;
#endif
				const Cell next = swipe_cell_update(it.diag(), match_scores, cbs_buf(i), extend_penalty, open_penalty, hgap, vgap, col_best, it.trace_mask(), row_counter, IdMask(p.query[i], target_seq));

				/*std::cout << "j=" << j << " i=" << i << " score=" << ScoreTraits<_sv>::int_score(extract_channel(next, 0)) <<
					" q=" << value_traits.alphabet[query[i]] << " t=" << value_traits.alphabet[extract_channel(target_seq, 0)]
					<< extract_stats(next, 0) << std::endl;*/

				it.set_hgap(hgap);
				it.set_score(next);
				++it;
			}
#ifdef STRICT_BAND
		}
#endif

		Score col_best_[CHANNELS], i_max[CHANNELS];
		store_sv(col_best, col_best_);
		row_counter.store(i_max);
		for (int i = 0; i < targets.active.size();) {
			int channel = targets.active[i];
			if (!targets.inc(channel))
				targets.active.erase(i);
			else
				++i;
			if (col_best_[channel] > best[channel]) {
				best[channel] = col_best_[channel];
				max_col[channel] = j;
				max_band_row[channel] = ScoreTraits<_sv>::int_score(i_max[channel]);
				stats[channel] = extract_stats(dp[max_band_row[channel]], channel);
				//std::cout << "stats[" << channel << "]=" << stats[channel] << std::endl;
			}
		}
		++i0;
		++i1;
		++j;
	}

	list<Hsp> out;
	TaskTimer timer;
	for (int i = 0; i < targets.n_targets; ++i) {
		if (best[i] < ScoreTraits<_sv>::max_score() && !overflow_stats<_sv>(stats[i])) {
			int score = ScoreTraits<_sv>::int_score(best[i]);
			if (!subject_begin[i].adjusted_matrix())
				score *= config.cbs_matrix_scale;
			const double evalue = score_matrix.evalue(score, qlen, subject_begin[i].true_target_len);
			if (score > 0 && score_matrix.report_cutoff(score, evalue)) {
				out.push_back(traceback<_sv>(composition_bias, dp, subject_begin[i], d_begin[i], best[i], evalue, max_col[i], i, i0 - j, i1 - j, max_band_row[i], stats[i], p));
			}
		}
		else
			overflow.push_back(subject_begin[i]);
	}
	p.stat.inc(Statistics::TIME_TRACEBACK, timer.microseconds());
	return out;
}

}}}