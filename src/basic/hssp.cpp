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

interval untranslate_range(const interval &r, unsigned frame, size_t len)
{
	signed f = frame <= 2 ? frame + 1 : 2 - frame;
	interval s;
	if (f > 0) {
		s.begin_ = (f - 1) + 3 * r.begin_;
		s.end_ = s.begin_ + 3 * r.length();
	}
	else {
		s.end_ = (int)len + f - 3 * r.begin_ + 1;
		s.begin_ = s.end_ - 3 * r.length();
	}
	return s;
}

bool Hsp_data::pass_through(const Diagonal_segment &d) const
{
	if (intersect(d.query_range(), query_range).length() != d.len
		|| intersect(d.subject_range(), subject_range).length() != d.len)
		return false;

	Iterator it = begin();
	const unsigned subject_end = d.j + d.len;
	const unsigned diag = d.diag();
	while (it.good()) {
		if ((int)it.subject_pos >= d.j) {
			if (it.subject_pos >= subject_end)
				return true;
			if (it.query_pos - it.subject_pos != diag)
				return false;
		}
		++it;
	}
	return true;
}

std::pair<int, int> Hsp_data::diagonal_bounds() const
{
	int d0 = std::numeric_limits<int>::max(), d1 = std::numeric_limits<int>::min();
	for (Iterator it = begin(); it.good(); ++it) {
		const int d = (int)it.query_pos - (int)it.subject_pos;
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

void Hsp_data::set_source_range(unsigned frame, unsigned dna_len)
{
	this->frame = frame;
	if (!align_mode.query_translated)
		query_source_range = query_range;
	else
		query_source_range = untranslate_range(query_range, frame, dna_len);
}

Hsp_context& Hsp_context::parse()
{
	hsp_.length = hsp_.identities = hsp_.mismatches = hsp_.gap_openings = hsp_.positives = hsp_.gaps = 0;
	unsigned d = 0;
	Iterator i = begin();

	for (; i.good(); ++i) {
		++hsp_.length;
		assert(i.query_pos < query.length());
		if (i.query_pos >= query.length())
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
		}
	}

	hsp_.query_range.end_ = i.query_pos;
	hsp_.subject_range.end_ = i.subject_pos;

	return *this;
}

Hsp_context& Hsp_context::set_query_source_range(unsigned oriented_query_begin)
{
	if(align_mode.query_translated)
		hsp_.set_source_range(oriented_query_begin);
	else
		hsp_.query_source_range = hsp_.query_range;
	return *this;
}