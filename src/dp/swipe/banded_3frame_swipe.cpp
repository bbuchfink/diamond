/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include <algorithm>
#include <thread>
#include "../dp.h"
#include "swipe_matrix.h"
#include "swipe.h"
#include "target_iterator.h"
#include "../../util/thread.h"

using namespace std;

template<typename _score> thread_local std::vector<score_vector<_score>> BandedSwipeMatrix<_score>::hgap_;
template<typename _score> thread_local std::vector<score_vector<_score>> BandedSwipeMatrix<_score>::score_;
template<typename _score> thread_local std::vector<score_vector<_score>> BandedSwipeTracebackMatrix<_score>::hgap_;
template<typename _score> thread_local std::vector<score_vector<_score>> BandedSwipeTracebackMatrix<_score>::score_;
template<typename _sv> thread_local std::vector<_sv> Banded3FrameSwipeMatrix<_sv>::hgap_;
template<typename _sv> thread_local std::vector<_sv> Banded3FrameSwipeMatrix<_sv>::score_;
template<typename _sv> thread_local std::vector<_sv> Banded3FrameSwipeTracebackMatrix<_sv>::hgap_;
template<typename _sv> thread_local std::vector<_sv> Banded3FrameSwipeTracebackMatrix<_sv>::score_;


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

template<typename _sv>
void traceback(sequence *query, Strand strand, int dna_len, const Banded3FrameSwipeTracebackMatrix<_sv> &dp, DpTarget &target, typename ScoreTraits<_sv>::Score max_score, int max_col, int channel, int i0, int i1, bool parallel)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	const int j0 = i1 - (target.d_end - 1), d1 = target.d_end, d0 = target.d_begin;
	typename Banded3FrameSwipeTracebackMatrix<_sv>::TracebackIterator it(dp.traceback(max_col + 1, i0 + max_col, j0 + max_col, dna_len, channel, max_score));
	if (parallel)
		target.tmp = new Hsp();
	else
		target.out->emplace_back();
		
	Hsp &out = parallel ? *target.tmp : target.out->back();
	out.score = target.score = ScoreTraits<_sv>::int_score(max_score);
	out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));

	out.set_end(it.i + 1, it.j + 1, Frame(strand, it.frame), dna_len);

	while (it.score() > ScoreTraits<_sv>::zero_score()) {
		const Letter q = query[it.frame][it.i], s = target.seq[it.j];
		const Score m = score_matrix(q, s), score = it.score();
		if (score == it.sm3() + m) {
			out.push_match(q, s, m > (Score)0);
			it.walk_diagonal();
		}
		else if (score == it.sm4() + m - score_matrix.frame_shift()) {
			out.push_match(q, s, m > (Score)0);
			out.transcript.push_back(op_frameshift_forward);
			it.walk_forward_shift();
		}
		else if (score == it.sm2() + m - score_matrix.frame_shift()) {
			out.push_match(q, s, m > (Score)0);
			out.transcript.push_back(op_frameshift_reverse);
			it.walk_reverse_shift();
		}
		else {
			const pair<Edit_operation, int> g(it.walk_gap(d0, d1));
			out.push_gap(g.first, g.second, &target.seq[it.j + g.second]);
		}
	}

	out.set_begin(it.i + 1, it.j + 1, Frame(strand, it.frame), dna_len);
	out.transcript.reverse();
	out.transcript.push_terminator();
}

template<typename _sv>
void traceback(sequence *query, Strand strand, int dna_len, const Banded3FrameSwipeMatrix<_sv> &dp, DpTarget &target, typename ScoreTraits<_sv>::Score max_score, int max_col, int channel, int i0, int i1, bool parallel)
{
	if (parallel)
		target.tmp = new Hsp();
	else
		target.out->emplace_back();
	Hsp &out = parallel ? *target.tmp : target.out->back();
	
	const int j0 = i1 - (target.d_end - 1);
	out.score = target.score = ScoreTraits<_sv>::int_score(max_score);
	out.query_range.end_ = std::min(i0 + max_col + (int)dp.band() / 3 / 2, (int)query[0].length());
	out.query_range.begin_ = std::max(out.query_range.end_ - (j0 + max_col), 0);
	out.frame = strand == FORWARD ? 0 : 3;
	out.query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(out.query_range.begin_, Frame(out.frame)), TranslatedPosition(out.query_range.end_, Frame(out.frame)), dna_len);
}

template<typename _sv, typename _traceback>
void banded_3frame_swipe(const TranslatedSequence &query, Strand strand, vector<DpTarget>::iterator subject_begin, vector<DpTarget>::iterator subject_end, DpStat &stat, bool parallel)
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
	
	for (int i = 0; i < targets.n_targets; ++i) {
		if (best[i] < ScoreTraits<_sv>::max_score()) {
			subject_begin[i].overflow = false;
			traceback<_sv>(q, strand, (int)query.source().length(), dp, subject_begin[i], best[i], max_col[i], i, i0 - j, i1 - j, parallel);
		}
		else
			subject_begin[i].overflow = true;
	}
}

void banded_3frame_swipe_worker(vector<DpTarget>::iterator begin,
	vector<DpTarget>::iterator end,
	Atomic<size_t> *next,
	bool score_only,
	const TranslatedSequence *query,
	Strand strand)
{
	DpStat stat;
	size_t pos;
	while (begin + (pos = next->post_add(config.swipe_chunk_size)) < end) {
		const auto e = min(begin + pos + config.swipe_chunk_size, end);
		for (vector<DpTarget>::iterator i = begin + pos; i < e; i += 8) {
			if (score_only || config.disable_traceback)
				banded_3frame_swipe<score_vector<int16_t>, ScoreOnly>(*query, strand, i, min(i + 8, end), stat, true);
			else
				banded_3frame_swipe<score_vector<int16_t>, Traceback>(*query, strand, i, min(i + 8, end), stat, true);
		}
	}
}

void banded_3frame_swipe(const TranslatedSequence &query, Strand strand, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end, DpStat &stat, bool score_only, bool parallel)
{
#ifdef __SSE2__
	task_timer timer("Banded 3frame swipe (sort)", parallel ? 3 : UINT_MAX);
	std::stable_sort(target_begin, target_end);
	if (parallel) {
		timer.go("Banded 3frame swipe (run)");
		vector<thread> threads;
		Atomic<size_t> next(0);
		for (size_t i = 0; i < config.threads_; ++i)
			threads.emplace_back(banded_3frame_swipe_worker,
				target_begin,
				target_end,
				&next,
				score_only,
				&query,
				strand);
		for (auto &t : threads)
			t.join();
		timer.go("Banded 3frame swipe (merge)");
		for (auto i = target_begin; i < target_end; ++i) {
			i->out->push_back(*i->tmp);
			delete i->tmp;
		}
	}
	else {
		for (vector<DpTarget>::iterator i = target_begin; i < target_end; i += 8)
			if (score_only || config.disable_traceback)
				banded_3frame_swipe<score_vector<int16_t>, ScoreOnly>(query, strand, i, std::min(i + 8, target_end), stat, false);
			else
				banded_3frame_swipe<score_vector<int16_t>, Traceback>(query, strand, i, std::min(i + 8, target_end), stat, false);
	}

	for (vector<DpTarget>::iterator i = target_begin; i < target_end; ++i)
		if (i->overflow) {
			if (score_only || config.disable_traceback)
				banded_3frame_swipe<int32_t, ScoreOnly>(query, strand, i, i + 1, stat, false);
			else
				banded_3frame_swipe<int32_t, Traceback>(query, strand, i, i + 1, stat, false);
		}
#else
	for (vector<DpTarget>::iterator i = target_begin; i < target_end; ++i)
		if (score_only || config.disable_traceback)
			banded_3frame_swipe<int32_t, ScoreOnly>(query, strand, i, i + 1, stat, false);
		else
			banded_3frame_swipe<int32_t, Traceback>(query, strand, i, i + 1, stat, false);
#endif
}