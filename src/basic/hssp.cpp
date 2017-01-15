/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include "../align/align.h"

bool Hsp_data::pass_through(const Diagonal_segment &d) const
{
	if (intersect(d.query_range(), query_range).length() != (size_t)d.len
		|| intersect(d.subject_range(), subject_range).length() != (size_t)d.len)
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
	else {
		signed f = frame <= 2 ? frame + 1 : 2 - frame;
		if (f > 0) {
			query_source_range.begin_ = (f - 1) + 3 * query_range.begin_;
			query_source_range.end_ = query_source_range.begin_ + 3 * query_range.length();
			//query_end_dna = (f-1) + 3 * (l.query_begin_+l.query_len_-1) + 3;
		}
		else {
			query_source_range.end_ = dna_len + f - 3 * query_range.begin_ + 1;
			//query_end_dna = dna_len + (f + 1) - 3 * (l.query_begin_+l.query_len_-1) - 2;
			query_source_range.begin_ = query_source_range.end_ - 3 * query_range.length();
		}
	}
}

Hsp_context& Hsp_context::parse()
{
	hsp_.length = hsp_.identities = hsp_.mismatches = hsp_.gap_openings = hsp_.positives = hsp_.gaps = 0;
	unsigned d = 0;
	Iterator i = begin();

	for (; i.good(); ++i) {
		++hsp_.length;
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