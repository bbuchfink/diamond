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

#include "../../data/reference.h"

template<typename _score> TLS_PTR vector<score_vector<_score> >* Banded3FrameSwipeTracebackMatrix<_score>::hgap_ptr;
template<typename _score> TLS_PTR vector<score_vector<_score> >* Banded3FrameSwipeTracebackMatrix<_score>::score_ptr;

#ifdef __SSE2__

template<typename _score, bool _with_transcript>
void traceback(sequence *query, Strand strand, int dna_len, const Banded3FrameSwipeTracebackMatrix<_score> &dp, DpTarget &target, _score max_score, int max_col, int channel, int i0, int i1)
{
	const int j0 = i1 - (target.d_end - 1), d1 = target.d_end, d0 = target.d_begin;
	typename Banded3FrameSwipeTracebackMatrix<_score>::TracebackIterator it(dp.traceback(max_col + 1, i0 + max_col, j0 + max_col, dna_len, channel, max_score));
	target.out->push_back(Hsp());
	Hsp &out = target.out->back();
	out.score = target.score = (uint16_t)max_score ^ 0x8000;
	if(_with_transcript) out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));

	out.set_end(it.i + 1, it.j + 1, Frame(strand, it.frame), dna_len);

	while (it.score() > SHRT_MIN) {
		const Letter q = query[it.frame][it.i], s = target.seq[it.j];
		const int m = score_matrix(q, s);
		const _score score = it.score();
		if (score == it.sm3() + m) {
			if (_with_transcript) {
				if (q == s)
					out.transcript.push_back(op_match, 1u);
				else
					out.transcript.push_back(op_substitution, s);
			}
			it.walk_diagonal();
		}
		else if (score == it.sm4() + m - score_matrix.frame_shift()) {
			if (_with_transcript) {
				if (q == s)
					out.transcript.push_back(op_match, 1u);
				else
					out.transcript.push_back(op_substitution, s);
				out.transcript.push_back(op_frameshift_forward);
			}
			it.walk_forward_shift();
		}
		else if (score == it.sm2() + m - score_matrix.frame_shift()) {
			if (_with_transcript) {
				if (q == s)
					out.transcript.push_back(op_match, 1u);
				else
					out.transcript.push_back(op_substitution, s);
				out.transcript.push_back(op_frameshift_reverse);
			}
			it.walk_reverse_shift();
		}
		else {
			int l = it.walk_gap(d0, d1);
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
		++out.length;
	}

	out.set_begin(it.i + 1, it.j + 1, Frame(strand, it.frame), dna_len);
	if (_with_transcript) {
		out.transcript.reverse();
		out.transcript.push_terminator();
	}
}

template<typename _score>
void banded_3frame_swipe(const TranslatedSequence &query, Strand strand, vector<DpTarget>::iterator subject_begin, vector<DpTarget>::iterator subject_end)
{
	typedef score_vector<_score> sv;

	assert(subject_end - subject_begin <= score_traits<_score>::channels);
	sequence q[3];
	query.get_strand(strand, q);
	const int qlen = (int)q[0].length(), qlen2 = (int)q[1].length(), qlen3 = (int)q[2].length();

	int band = 0;
	for (vector<DpTarget>::const_iterator j = subject_begin; j < subject_end; ++j) {
		band = std::max(band, j->d_end - j->d_begin);
	}

	int i0 = INT_MAX, i1 = INT_MAX;
	for (vector<DpTarget>::iterator j = subject_begin; j < subject_end; ++j) {
		if (j->d_end - j->d_begin < band) {
			const int top_max = j->d_begin - (-(int)(j->seq.length() - 1)), bottom_max = qlen - j->d_end;
			int diff = band - (j->d_end - j->d_begin);
			int d = std::min(std::max(diff / 2, diff - bottom_max), top_max);
			j->d_begin -= d;
			diff -= d;
			if (diff > bottom_max)
				throw std::runtime_error("xxx");
			j->d_end += std::min(diff, bottom_max);
		}

		int i2 = std::max(j->d_end - 1, 0);
		i1 = std::min(i1, i2);
		i0 = std::min(i0, i2 + 1 - band);
	}

	TargetIterator<score_traits<_score>::channels> targets(subject_begin, subject_end, i1, qlen);
	Banded3FrameSwipeTracebackMatrix<_score> dp(band * 3, targets.cols);

	const sv open_penalty(score_matrix.gap_open() + score_matrix.gap_extend()),
		extend_penalty(score_matrix.gap_extend()),
		frameshift_penalty(score_matrix.frame_shift());
	
	SwipeProfile<_score> profile;
	int16_t best[8];
	int max_col[8];
	for (int i = 0; i < 8; ++i)
		best[i] = SHRT_MIN;

	int j = 0;
	while (targets.active.size() > 0) {
		const int i0_ = std::max(i0, 0), i1_ = std::min(i1, qlen - 1);
		if (i0_ > i1_)
			break;
		typename Banded3FrameSwipeTracebackMatrix<_score>::ColumnIterator it(dp.begin((i0_ - i0) * 3, j));
		if (i0_ - i0 > 0)
			it.set_zero();
		sv vgap0, vgap1, vgap2, hgap;
		vgap0.zero();
		vgap1.zero();
		vgap2.zero();
		sv col_best;
		col_best.zero();

		profile.set(targets.get());
		for (int i = i0_; i <= i1_; ++i) {
			hgap = it.hgap();
			sv next = cell_update<_score>(it.sm3, it.sm4, it.sm2, profile.get(q[0][i]), extend_penalty, open_penalty, frameshift_penalty, hgap, vgap0, col_best);
			it.set_hgap(hgap);
			it.set_score(next);
			++it;

			if (i >= qlen2)
				break;
			hgap = it.hgap();
			next = cell_update<_score>(it.sm3, it.sm4, it.sm2, profile.get(q[1][i]), extend_penalty, open_penalty, frameshift_penalty, hgap, vgap1, col_best);
			it.set_hgap(hgap);
			it.set_score(next);
			++it;

			if (i >= qlen3)
				break;
			hgap = it.hgap();
			next = cell_update<_score>(it.sm3, it.sm4, it.sm2, profile.get(q[2][i]), extend_penalty, open_penalty, frameshift_penalty, hgap, vgap2, col_best);
			it.set_hgap(hgap);
			it.set_score(next);
			++it;
		}

		cells += targets.active.size()*(i1_ - i0_ + 1) * 3;
		for (int i = 0; i < targets.active.size();) {
			int channel = targets.active[i];
			if (!targets.inc(channel))
				targets.active.erase(i);
			else
				++i;
			if (col_best[channel] > best[channel]) {
				best[channel] = col_best[channel];
				max_col[channel] = j;
			}
		}
		++i0;
		++i1;
		++j;
	}
	for (int i = 0; i < targets.n_targets; ++i)
		traceback<_score, true>(q, strand, (int)query.source().length(), dp, subject_begin[i], best[i], max_col[i], i, i0 - j, i1 - j);
}

#endif

void banded_3frame_swipe(const TranslatedSequence &query, Strand strand, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end)
{
	/*cout << "==========" << endl;
	for (vector<DpTarget>::iterator i = target_begin; i < target_end; ++i)
		cout << ref_ids::get()[i->subject_id].c_str() << endl;*/
#ifdef __SSE2__
	banded_3frame_swipe<int16_t>(query, strand, target_begin, target_end);
#endif
}