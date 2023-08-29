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

#include <vector>
#include "swipe.h"
#include "../../basic/sequence.h"
#include "target_iterator.h"
#include "../../util/data_structures/mem_buffer.h"

using std::vector;
using std::list;
using std::max;
using namespace DISPATCH_ARCH;

namespace DP { namespace Swipe {
namespace DISPATCH_ARCH {

template<typename _sv, typename Cell, typename _cbs, typename StatType>
Hsp traceback(_cbs bias_correction, const Matrix<Cell>& dp, const DpTarget& target, typename ScoreTraits<_sv>::Score max_score, double evalue, int max_col, int max_i, int max_j, int channel, const StatType &stats, Params& p)
{
	Hsp out(false);
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<_sv>::int_score(max_score) * config.cbs_matrix_scale;
	out.evalue = evalue;
	out.bit_score = score_matrix.bitscore(out.score);
	out.corrected_bit_score = score_matrix.bitscore_corrected(out.score, p.query.length(), target.true_target_len);
	out.frame = p.frame.index();
	if (target.carry_over.i1 == 0) {
		out.query_range.end_ = max_i + 1;
		out.subject_range.end_ = max_j + 1;
	}
	else {
		out.query_range.end_ = target.carry_over.i1;
		out.subject_range.end_ = target.carry_over.j1;
		out.identities = target.carry_over.ident;
		out.length = target.carry_over.len;
		out.query_range.begin_ = (int)p.query.length() - 1 - max_i;
		out.subject_range.begin_ = (int)target.seq.length() - 1 - max_j;
		try {
			out.approx_id = out.approx_id_percent(Sequence(p.query.reverse()), Sequence(target.seq.reverse()));
		}
		catch (std::out_of_range&) {
			throw std::runtime_error(std::string("Out_of_range query=") + std::string(p.query_id) + " target=" + target.seq.to_string());
		}
	}
	out.target_seq = target.seq;
	out.matrix = target.matrix;
	assign_stats(out, stats);
	out.query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(out.query_range.begin_, p.frame), TranslatedPosition(out.query_range.end_, p.frame), p.query_source_len);
	return out;
}

template<typename _sv, typename _cbs>
Hsp traceback(_cbs bias_correction, const TracebackVectorMatrix<_sv> &dp, const DpTarget &target, typename ScoreTraits<_sv>::Score max_score, double evalue, int max_col, int max_i, int max_j, int channel, Void, Params& p)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	typedef typename ScoreTraits<_sv>::TraceMask TraceMask;
	const auto channel_mask = TraceMask::vmask(channel) | TraceMask::hmask(channel);
	typename TracebackVectorMatrix<_sv>::TracebackIterator it(dp.traceback(max_col, max_i, max_j, channel));
	Hsp out(true);
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<_sv>::int_score(max_score);
	out.evalue = evalue;
	out.bit_score = score_matrix.bitscore(out.score);
	out.corrected_bit_score = score_matrix.bitscore_corrected(out.score, p.query.length(), target.true_target_len);
	out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));

	out.frame = p.frame.index();
	out.query_range.end_ = it.i + 1;
	out.subject_range.end_ = it.j + 1;
	const int end_score = out.score;
	int score = 0;
	const bool adjusted_matrix = target.adjusted_matrix();
	if (!adjusted_matrix)
		out.score *= config.cbs_matrix_scale;
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
			score -= score_matrix.gap_open() + g.second * score_matrix.gap_extend();
		}
	}

	if (score != end_score)
		throw std::runtime_error("Traceback error. " + p.query.to_string());

	out.query_range.begin_ = it.i + 1;
	out.subject_range.begin_ = it.j + 1;
	out.transcript.reverse();
	out.transcript.push_terminator();
	out.query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(out.query_range.begin_, p.frame), TranslatedPosition(out.query_range.end_, p.frame), p.query_source_len);
	out.approx_id = out.approx_id_percent(p.query, target.seq);
	return out;
}

template<typename _sv, typename _cbs, typename It, typename Cfg>
list<Hsp> swipe(const It target_begin, const It target_end, std::atomic<BlockId>* const next, _cbs composition_bias, vector<DpTarget>& overflow, Params& p)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	using Cell = typename Cfg::Cell;
	typedef typename SelectMatrix<Cell, Cfg::traceback>::Type Matrix;
	using RowCounter = typename Cfg::RowCounter;
	using IdMask = typename Cfg::IdMask;
	using StatType = decltype(extract_stats<_sv>(Cell(), int()));
	constexpr int CHANNELS = ScoreTraits<_sv>::CHANNELS;

	int max_col[CHANNELS], max_i[CHANNELS], max_j[CHANNELS];
	const int qlen = (int)p.query.length();

	if (qlen > RowCounter::MAX_LEN)
		throw std::runtime_error("Query length exceeds row counter maximum.");

	if (config.cbs_matrix_scale != 1)
		throw std::runtime_error("Matrix scale != 1.0 not supported.");

	const _sv open_penalty(static_cast<Score>(score_matrix.gap_open() + score_matrix.gap_extend())),
		extend_penalty(static_cast<Score>(score_matrix.gap_extend()));
	//_sv best = _sv();
	Score best[CHANNELS];
	StatType hsp_stats[CHANNELS];
	std::fill(best, best + CHANNELS, ScoreTraits<_sv>::zero_score());
	SwipeProfile<_sv> profile;
	std::array<const int8_t*, 32> target_scores;
	AsyncTargetBuffer<Score, It> targets(target_begin, target_end, next);
	Matrix dp(qlen, targets.max_len());
	CBSBuffer<_sv, _cbs> cbs_buf(composition_bias, qlen, 0);
	list<Hsp> out;
	int col = 0;
	
	while (targets.active.size() > 0) {
		typename Matrix::ColumnIterator it(dp.begin(col));
		RowCounter row_counter(0);
		Cell vgap, hgap, last;
		_sv col_best;
		vgap = hgap = last = col_best = _sv();

		const auto target_seq_vector = targets.seq_vector();
		const _sv target_seq(target_seq_vector);
		if (targets.cbs_mask() != 0) {
			if (targets.custom_matrix_16bit)
				profile.set(targets.get32().data());
			else
				profile.set(targets.get(target_scores.data()));
		}
		else {
#ifdef __SSSE3__
			profile.set(target_seq_vector);
#else
			profile.set(targets.get(target_scores.data()));
#endif
		}

#ifdef DP_STAT
		p.stat.inc(Statistics::GROSS_DP_CELLS, uint64_t(qlen) * CHANNELS);
#endif
		for (int i = 0; i < qlen; ++i) {
			hgap = it.hgap();
			const Cell next = swipe_cell_update(it.diag(), profile.get(p.query[i]), cbs_buf(i), extend_penalty, open_penalty, hgap, vgap, col_best, it.trace_mask(), row_counter, IdMask(p.query[i], target_seq));

			/*/std::cout << "j=" << targets.pos[0] << " i=" << i << " score=" << ScoreTraits<_sv>::int_score(extract_channel(next, 0)) <<
				" q=" << value_traits.alphabet[p.query[i]] << " t=" << value_traits.alphabet[extract_channel(target_seq, 0)]
				<< extract_stats(next, 0) << " j'=" << targets.dp_targets[0].seq.length() - targets.pos[0] << " i'=" << p.query.length() - i << std::endl;*/

			it.set_hgap(hgap);
			it.set_score(last);
			last = next;
			++it;
		}
		it.set_score(last);
		//best = max(best, col_best);
		
		Score col_best_[CHANNELS], i_max[CHANNELS];
		store_sv(col_best, col_best_);
		row_counter.store(i_max);
		for (int i = 0; i < targets.active.size();) {
			int c = targets.active[i];
			if (col_best_[c] > best[c]) {
				best[c] = col_best_[c];
				max_col[c] = col;
				max_i[c] = ScoreTraits<_sv>::int_score(i_max[c]);
				max_j[c] = targets.pos[c];
				hsp_stats[c] = extract_stats(dp[max_i[c]], c);
				//std::cout << "stats[" << c << "]=" << hsp_stats[c] << " j=" << targets.pos[0] << " j'=" << targets.dp_targets[0].seq.length() - targets.pos[0] << " score=" << ScoreTraits<_sv>::int_score(best[c]) << std::endl;
			}
			bool reinit = false;
			if (col_best_[c] == ScoreTraits<_sv>::max_score()) {
				overflow.push_back(targets.dp_targets[c]);
				reinit = true;
			} else if (!targets.inc(c)) {
				if (overflow_stats<_sv>(hsp_stats[c]))
					overflow.push_back(targets.dp_targets[c]);
				else {
					const int s = ScoreTraits<_sv>::int_score(best[c]) * config.cbs_matrix_scale;
					const double evalue = score_matrix.evalue(s, qlen, (unsigned)targets.dp_targets[c].true_target_len);
					if (s > 0 && score_matrix.report_cutoff(s, evalue))
						out.push_back(traceback<_sv>(composition_bias, dp, targets.dp_targets[c], best[c], evalue, max_col[c], max_i[c], max_j[c], c, hsp_stats[c], p));
				}
				reinit = true;
			}
			if (reinit) {
				if (targets.init_target(i, c)) {
					dp.set_zero(c);
					//set_channel(best, c, ScoreTraits<_sv>::zero_score());
					best[c] = ScoreTraits<_sv>::zero_score();
				}
				else
					continue;
			}
			++i;
		}
		col = (col + 1) % dp.cols();
	}

	return out;
}

}}}