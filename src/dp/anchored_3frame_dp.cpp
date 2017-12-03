/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#include "dp.h"
#include "last/GappedXdropAligner.hh"

void anchored_3frame_dp(const TranslatedSequence &query, sequence &subject, const DiagonalSegment &anchor, Hsp_data &out, int gap_open, int gap_extend, int frame_shift)
{
	static TLS_PTR cbrc::GappedXdropAligner* al_ptr;
	cbrc::GappedXdropAligner al(TLS::get(al_ptr));

	Strand strand = anchor.i.frame.strand;

	const unsigned char* q[3];
	out.clear();

	// Extend left

	q[0] = (unsigned char*)&query[anchor.i];
	TranslatedPosition forward(anchor.i);
	forward.shift_forward();
	q[2] = (unsigned char*)&query[forward];
	TranslatedPosition back(anchor.i);
	back.shift_back();
	q[1] = (unsigned char*)&query[back];

	out.score = al.align3((unsigned char*)&subject[anchor.j],
		q[0],
		q[1],
		q[2],
		false,
		(const cbrc::ScoreMatrixRow*)score_matrix.matrix32(),
		gap_open,
		gap_extend,
		cbrc::INF,
		frame_shift,
		score_matrix.rawscore(config.gapped_xdrop),
		score_matrix.high_score());

	size_t end1, end2, length;
	DiagonalSegment last;

	while (al.getNextChunk3(end1,
		end2,
		length,
		gap_open,
		gap_extend,
		cbrc::INF,
		frame_shift)) {
		//cout << "end1=" << end1 << " end2=" << end2 << " length=" << length << endl;
		const DiagonalSegment d(TranslatedPosition(anchor.i.in_strand() - end2, strand), anchor.j - end1, (int)length);
		if (last.len > 0)
			out.splice(last, d, query, subject, false);
		else
			out.set_begin(d, query.source().length());
		out.push_back(d, query, subject, false);
		last = d;
	}

	// Splice with anchor
	//int splice_score = last.splice_score(anchor, gap_open, gap_extend, frame_shift);
	if (last.len == 0) {
		out.set_begin(anchor, query.source().length());
	}
	else
		out.splice(last, anchor, query, subject, false);

	// Push anchor
	out.push_back(anchor, query, subject, false);
	out.score += anchor.score;
	const size_t splice = out.transcript.raw_length();
	out.transcript.push_terminator();

	// Extend right
	
	q[0] = (unsigned char*)&query[anchor.query_end()];
	forward = anchor.query_end();
	forward.shift_forward();
	q[1] = (unsigned char*)&query[forward];
	back = anchor.query_end();
	back.shift_back();
	q[2] = (unsigned char*)&query[back];

	out.score += al.align3((unsigned char*)&subject[anchor.subject_end()],
		q[0],
		q[1],
		q[2],
		true,
		(const cbrc::ScoreMatrixRow*)score_matrix.matrix32(),
		gap_open,
		gap_extend,
		cbrc::INF,
		frame_shift,
		score_matrix.rawscore(config.gapped_xdrop),
		score_matrix.high_score());

	last = DiagonalSegment();

	while (al.getNextChunk3(end1,
		end2,
		length,
		gap_open,
		gap_extend,
		cbrc::INF,
		frame_shift)) {
		//cout << "end1=" << end1 << " end2=" << end2 << " length=" << length << endl;
		const DiagonalSegment d(TranslatedPosition(anchor.query_end().in_strand() + end2 - length * 3, strand), anchor.subject_end() + end1 - length, (int)length);
		if (last.len > 0)
			out.splice(d, last, query, subject, true);
		else
			out.set_end(d, query.source().length());
		out.push_back(d, query, subject, true);
		last = d;
	}

	// Splice with anchor
	//slice_score = anchor.splice_score(last, gap_open, gap_extend, frame_shift);
	if (last.len == 0) {
		out.set_end(anchor, query.source().length());
	}
	else
		out.splice(anchor, last, query, subject, true);

	out.transcript.reverse(splice);
	out.transcript.push_terminator();
}