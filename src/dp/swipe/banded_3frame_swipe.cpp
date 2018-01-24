/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "../dp.h"
#include "swipe_matrix.h"
#include "swipe.h"
#include "target_iterator.h"

// #include "../../data/reference.h"

template<typename _sv> TLS_PTR vector<_sv>* Banded3FrameSwipeTracebackMatrix<_sv>::hgap_ptr;
template<typename _sv> TLS_PTR vector<_sv>* Banded3FrameSwipeTracebackMatrix<_sv>::score_ptr;
template<typename _sv> TLS_PTR vector<_sv>* Banded3FrameSwipeMatrix<_sv>::hgap_ptr;
template<typename _sv> TLS_PTR vector<_sv>* Banded3FrameSwipeMatrix<_sv>::score_ptr;

struct Traceback {};
struct ScoreOnly {};

template<typename _sv, typename _traceback>
struct Banded3FrameSwipeMatrixRef
{
};

template<typename _sv>
struct Banded3FrameSwipeMatrixRef<_sv, Traceback>
{
	typedef Banded3FrameSwipeTracebackMatrix<_sv> type;
};

template<typename _sv>
struct Banded3FrameSwipeMatrixRef<_sv, ScoreOnly>
{
	typedef Banded3FrameSwipeMatrix<_sv> type;
};

template<typename _sv, bool _with_transcript>
void traceback(sequence *query, Strand strand, int dna_len, const Banded3FrameSwipeTracebackMatrix<_sv> &dp, DpTarget &target, typename ScoreTraits<_sv>::Score max_score, int max_col, int channel, int i0, int i1)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	const int j0 = i1 - (target.d_end - 1), d1 = target.d_end, d0 = target.d_begin;
	typename Banded3FrameSwipeTracebackMatrix<_sv>::TracebackIterator it(dp.traceback(max_col + 1, i0 + max_col, j0 + max_col, dna_len, channel, max_score));
	target.out->push_back(Hsp());
	Hsp &out = target.out->back();
	out.score = target.score = ScoreTraits<_sv>::int_score(max_score);
	if(_with_transcript) out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));

	out.set_end(it.i + 1, it.j + 1, Frame(strand, it.frame), dna_len);

	while (it.score() > ScoreTraits<_sv>::zero_score()) {
		const Letter q = query[it.frame][it.i], s = target.seq[it.j];
		const Score m = score_matrix(q, s), score = it.score();
		if (score == it.sm3() + m) {
			if (_with_transcript) out.push_match(q, s, m > (Score)0);
			it.walk_diagonal();
		}
		else if (score == it.sm4() + m - score_matrix.frame_shift()) {
			if (_with_transcript) {
				out.push_match(q, s, m > (Score)0);
				out.transcript.push_back(op_frameshift_forward);
			}
			it.walk_forward_shift();
		}
		else if (score == it.sm2() + m - score_matrix.frame_shift()) {
			if (_with_transcript) {
				out.push_match(q, s, m > (Score)0);
				out.transcript.push_back(op_frameshift_reverse);
			}
			it.walk_reverse_shift();
		}
		else {
			int l = it.walk_gap(d0, d1);
			++out.gap_openings;
			out.length += abs(l);
			out.gaps += abs(l);
			if (_with_transcript) {
				if (l > 0) {
					out.transcript.push_back(op_insertion, (unsigned)l);
				}
				else {
					for (int j = it.j - l; j > it.j; --j)
						out.transcript.push_back(op_deletion, target.seq[j]);
				}
			}
		}
	}

	out.set_begin(it.i + 1, it.j + 1, Frame(strand, it.frame), dna_len);
	if (_with_transcript) {
		out.transcript.reverse();
		out.transcript.push_terminator();
	}
}

template<typename _sv, bool _with_transcript>
void traceback(sequence *query, Strand strand, int dna_len, const Banded3FrameSwipeMatrix<_sv> &dp, DpTarget &target, typename ScoreTraits<_sv>::Score max_score, int max_col, int channel, int i0, int i1)
{
	target.out->push_back(Hsp());
	Hsp &out = target.out->back();
	out.score = target.score = ScoreTraits<_sv>::int_score(max_score);
}

template<typename _sv, typename _traceback>
void banded_3frame_swipe(const TranslatedSequence &query, Strand strand, vector<DpTarget>::iterator subject_begin, vector<DpTarget>::iterator subject_end, DpStat &stat)
{
	typedef typename Banded3FrameSwipeMatrixRef<_sv, _traceback>::type Matrix;
	typedef typename ScoreTraits<_sv>::Score Score;

	assert(subject_end - subject_begin <= ScoreTraits<_sv>::CHANNELS);
	sequence q[3];
	query.get_strand(strand, q);
	const int qlen = (int)q[0].length(), qlen2 = (int)q[1].length(), qlen3 = (int)q[2].length();

	int band = 0;
	for (vector<DpTarget>::const_iterator j = subject_begin; j < subject_end; ++j)
		band = std::max(band, j->d_end - j->d_begin);

	int i0 = INT_MAX, i1 = INT_MAX;
	for (vector<DpTarget>::iterator j = subject_begin; j < subject_end; ++j) {
		/*if (j->d_end - j->d_begin < band) {
			const int top_max = j->d_begin - (-(int)(j->seq.length() - 1)), bottom_max = qlen - j->d_end;
			int diff = band - (j->d_end - j->d_begin);
			int d = std::min(std::max(diff / 2, diff - bottom_max), top_max);
			j->d_begin -= d;
			diff -= d;
			if (diff > bottom_max)
				throw std::runtime_error("");
			j->d_end += std::min(diff, bottom_max);
		}*/
		j->d_begin = j->d_end - band;
		int i2 = std::max(j->d_end - 1, 0);
		i1 = std::min(i1, i2);
		i0 = std::min(i0, i2 + 1 - band);
	}

	TargetIterator<ScoreTraits<_sv>::CHANNELS> targets(subject_begin, subject_end, i1, qlen);
	Matrix dp(band * 3, targets.cols);

	const _sv open_penalty(score_matrix.gap_open() + score_matrix.gap_extend()),
		extend_penalty(score_matrix.gap_extend()),
		frameshift_penalty(score_matrix.frame_shift());
	
	SwipeProfile<_sv> profile;
	Score best[ScoreTraits<_sv>::CHANNELS];
	int max_col[ScoreTraits<_sv>::CHANNELS];
	for (int i = 0; i < ScoreTraits<_sv>::CHANNELS; ++i)
		best[i] = ScoreTraits<_sv>::zero_score();

	int j = 0;
	while (targets.active.size() > 0) {
		const int i0_ = std::max(i0, 0), i1_ = std::min(i1, qlen - 1);
		if (i0_ > i1_)
			break;
		typename Matrix::ColumnIterator it(dp.begin((i0_ - i0) * 3, j));
		if (i0_ - i0 > 0)
			it.set_zero();
		_sv vgap0, vgap1, vgap2, hgap, col_best;
		vgap0 = vgap1 = vgap2 = col_best = ScoreTraits<_sv>::zero();

		profile.set(targets.get());
		for (int i = i0_; i <= i1_; ++i) {
			hgap = it.hgap();
			_sv next = cell_update<_sv>(it.sm3, it.sm4, it.sm2, profile.get(q[0][i]), extend_penalty, open_penalty, frameshift_penalty, hgap, vgap0, col_best);
			it.set_hgap(hgap);
			it.set_score(next);
			++it;

			if (i >= qlen2)
				break;
			hgap = it.hgap();
			next = cell_update<_sv>(it.sm3, it.sm4, it.sm2, profile.get(q[1][i]), extend_penalty, open_penalty, frameshift_penalty, hgap, vgap1, col_best);
			it.set_hgap(hgap);
			it.set_score(next);
			++it;

			if (i >= qlen3)
				break;
			hgap = it.hgap();
			next = cell_update<_sv>(it.sm3, it.sm4, it.sm2, profile.get(q[2][i]), extend_penalty, open_penalty, frameshift_penalty, hgap, vgap2, col_best);
			it.set_hgap(hgap);
			it.set_score(next);
			++it;
		}

#ifdef DP_STAT
		stat.net_cells += targets.live * (i1_ - i0_ + 1) * 3;
		stat.gross_cells += ScoreTraits<_sv>::CHANNELS * (i1_ - i0_ + 1) * 3;
#endif

		Score col_best_[ScoreTraits<_sv>::CHANNELS];
		store_sv(col_best, col_best_);
		for (int i = 0; i < targets.active.size();) {
			int channel = targets.active[i];
			if (!targets.inc(channel))
				targets.active.erase(i);
			else
				++i;
			if (col_best_[channel] > best[channel]) {
				best[channel] = col_best_[channel];
				max_col[channel] = j;
			}
		}
		++i0;
		++i1;
		++j;
	}
	for (int i = 0; i < targets.n_targets; ++i)
		if (best[i] < ScoreTraits<_sv>::max_score()) {
			subject_begin[i].overflow = false;
			traceback<_sv, true>(q, strand, (int)query.source().length(), dp, subject_begin[i], best[i], max_col[i], i, i0 - j, i1 - j);
		}
		else
			subject_begin[i].overflow = true;
}

void banded_3frame_swipe(const TranslatedSequence &query, Strand strand, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end, DpStat &stat, bool score_only)
{
	/*cout << "==========" << endl;
	for (vector<DpTarget>::iterator i = target_begin; i < target_end; ++i)
		cout << ref_ids::get()[i->subject_id].c_str() << endl;*/
	std::stable_sort(target_begin, target_end);
#ifdef __SSE2__
	for (vector<DpTarget>::iterator i = target_begin; i < target_end; i += 8)
		if (score_only || config.disable_traceback)
			banded_3frame_swipe<score_vector<int16_t>, ScoreOnly>(query, strand, i, std::min(i + 8, target_end), stat);
		else
			banded_3frame_swipe<score_vector<int16_t>, Traceback>(query, strand, i, std::min(i + 8, target_end), stat);
	for (vector<DpTarget>::iterator i = target_begin; i < target_end; ++i)
		if (i->overflow) {
			if (score_only || config.disable_traceback)
				banded_3frame_swipe<int32_t, ScoreOnly>(query, strand, i, i + 1, stat);
			else
				banded_3frame_swipe<int32_t, Traceback>(query, strand, i, i + 1, stat);
		}
#else
	for (vector<DpTarget>::iterator i = target_begin; i < target_end; i += 1)
		if (score_only || config.disable_traceback)
			banded_3frame_swipe<int32_t, ScoreOnly>(query, strand, i, i + 1, stat);
		else
			banded_3frame_swipe<int32_t, Traceback>(query, strand, i, i + 1, stat);
#endif
}