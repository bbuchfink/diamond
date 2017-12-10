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

#include "../align/align.h"

std::pair<int, int> Hsp_data::diagonal_bounds() const
{
	int d0 = std::numeric_limits<int>::max(), d1 = std::numeric_limits<int>::min();
	for (Iterator it = begin(); it.good(); ++it) {
		const int d = (int)it.query_pos.translated - (int)it.subject_pos;
		d0 = std::min(d0, d);
		d1 = std::max(d1, d);
	}
	return std::make_pair(d0, d1);
}

bool Hsp_data::is_weakly_enveloped(const Hsp_data &j) const
{
	static const double overlap_factor = 0.9;
	return score <= j.score
		&& subject_range.overlap_factor(j.subject_range) >= overlap_factor
		&& query_range.overlap_factor(j.query_range) >= overlap_factor;
}

void Hsp_data::merge(const Hsp_data &right, const Hsp_data &left, unsigned query_anchor, unsigned subject_anchor)
{
	length = right.length + left.length;
	gap_openings = right.gap_openings + left.gap_openings;
	identities = right.identities + left.identities;
	mismatches = right.mismatches + left.mismatches;
	score = right.score + left.score;
	gaps = right.gaps + left.gaps;
	gap_openings = right.gap_openings + left.gap_openings;
	positives = right.positives + left.positives;
	subject_range = interval(subject_anchor + 1 - left.subject_range.end_, subject_anchor + 1 + right.subject_range.end_);
	query_range = interval(query_anchor + 1 - left.query_range.end_, query_anchor + 1 + right.query_range.end_);
	transcript.data_.insert(transcript.data_.end(), left.transcript.data_.begin(), left.transcript.data_.end());
	transcript.data_.insert(transcript.data_.end(), right.transcript.data_.rbegin(), right.transcript.data_.rend());
	transcript.push_terminator();
}

Hsp_context& Hsp_context::parse()
{
	hsp_.length = hsp_.identities = hsp_.mismatches = hsp_.gap_openings = hsp_.positives = hsp_.gaps = 0;
	unsigned d = 0;
	Iterator i = begin();

	for (; i.good(); ++i) {
		++hsp_.length;
		assert(query.in_bounds(i.query_pos));
		if (!query.in_bounds(i.query_pos))
			throw std::runtime_error("Query sequence index out of bounds.");
		switch (i.op()) {
		case op_match:
			++hsp_.identities;
			++hsp_.positives;
			d = 0;
			break;
		case op_substitution:
			++hsp_.mismatches;
			if (i.score() > 0)
				++hsp_.positives;
			d = 0;
			break;
		case op_insertion:
		case op_deletion:
			if (d == 0)
				++hsp_.gap_openings;
			++d;
			++hsp_.gaps;
			break;
		default:
			;
		}
	}

	hsp_.query_range.end_ = i.query_pos.translated;
	hsp_.subject_range.end_ = i.subject_pos;
	hsp_.query_source_range = TranslatedPosition::absolute_interval(begin().query_pos, i.query_pos, query.source().length());

	return *this;
}

void Hsp_data::push_back(const DiagonalSegment &d, const TranslatedSequence &query, const sequence &subject, bool reversed)
{
	const sequence &q = query[d.i.frame];
	if (reversed) {
		for (int i = d.query_last().translated, j = d.subject_last(); j >= d.j; --j, --i) {
			const Letter ls = subject[j], lq = q[i];
			if (ls == lq) {
				transcript.push_back(op_match);
				++identities;
				++positives;
			}
			else {
				transcript.push_back(op_substitution, ls);
				++mismatches;
				if (score_matrix(ls, lq) > 0)
					++positives;
			}
			++length;
		}
	}
	else {
		for (int i = d.i, j = d.j; j < d.subject_end(); ++j, ++i) {
			const Letter ls = subject[j], lq = q[i];
			if (ls == lq) {
				transcript.push_back(op_match);
				++identities;
				++positives;
			}
			else {
				transcript.push_back(op_substitution, ls);
				++mismatches;
				if (score_matrix(ls, lq) > 0)
					++positives;
			}
			++length;
		}
	}
}

void Hsp_data::splice(const DiagonalSegment &a, const DiagonalSegment &b, const TranslatedSequence &query, const sequence &subject, bool reversed)
{
	TranslatedPosition i0 = a.query_last();
	int j0 = a.subject_last();
	const int fs = i0.frame_shift(b.i);
	if (fs == 1) {
		i0.shift_forward();
		transcript.push_back(op_frameshift_forward);
	}
	else if (fs == -1) {
		i0.shift_back();
		transcript.push_back(op_frameshift_reverse);
	}
	const int d0 = i0 - j0, d1 = b.diag();
	if (d1 > d0)
		transcript.push_back(op_insertion, unsigned(d1 - d0));
	else if (d1 < d0) {
		if (reversed)
			transcript.push_back(subject.subseq(j0 + 1, b.j), Reversed());
		else
			transcript.push_back(subject.subseq(j0 + 1, b.j));
	}
	const int shift = abs(d1 - d0);
	if (shift > 0) {
		length += shift;
		++gap_openings;
		gaps += shift;
	}
}

void Hsp_data::set_begin(const DiagonalSegment &d, int dna_len)
{
	subject_range.begin_ = d.j;
	query_range.begin_ = d.i;
	frame = d.i.frame.index();
	if (d.i.frame.strand == FORWARD)
		query_source_range.begin_ = d.i.absolute(dna_len);
	else
		query_source_range.end_ = d.i.absolute(dna_len) + 1;
}

void Hsp_data::set_end(const DiagonalSegment &d, int dna_len)
{
	subject_range.end_ = d.subject_end();
	query_range.end_ = d.query_end();
	if (d.i.frame.strand == FORWARD)
		query_source_range.end_ = d.query_end().absolute(dna_len);
	else
		query_source_range.begin_ = d.query_end().absolute(dna_len) + 1;
}

void Hsp_data::clear()
{
	score = frame = length = identities = mismatches = positives = gap_openings = gaps = sw_score = 0;
	transcript.clear();
}

bool Hsp_data::is_weakly_enveloped_by(list<Hsp_data>::const_iterator begin, list<Hsp_data>::const_iterator end, int cutoff) const
{
	for (list<Hsp_data>::const_iterator i = begin; i != end; ++i)
		if (partial_score(*i) < cutoff)
			return true;
	return false;
}