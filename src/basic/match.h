/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <list>
#include "sequence.h"
#include "value.h"
#include "packed_transcript.h"
#include "translated_position.h"
#include "util/geo/diagonal_segment.h"
#include "util/io/input_file.h"
#include "util/hsp/approx_hsp.h"

inline Interval normalized_range(unsigned pos, int len, Strand strand)
{
	return strand == FORWARD
			? Interval(pos, pos + len)
			: Interval(pos + 1 + len, pos + 1);
}

struct IntermediateRecord;
struct OutputFormat;

namespace Stats {
struct TargetMatrix;
struct Blastn_Score;
}

struct Hsp
{

	Hsp(const bool backtraced = false) :
		backtraced(backtraced),
		score(0),
		frame(0),
		length(0),
		identities(0),
		mismatches(0),
		positives(0),
		gap_openings(0),
		gaps(0),
		swipe_target(0),
		d_begin(0),
		d_end(0),
		evalue(DBL_MAX),
		bit_score(0.0),
		corrected_bit_score(0.0),
		approx_id(0.0),
#ifdef WITH_DNA
        mapping_quality(0),
        n_anchors(0),
#endif
		matrix(nullptr)
	{}

	Hsp(const bool backtraced, int score, int swipe_target = 0) :
		backtraced(backtraced),
		score(score),
		frame(0),
		length(0),
		identities(0),
		mismatches(0),
		positives(0),
		gap_openings(0),
		gaps(0),
		swipe_target(swipe_target),
		d_begin(0),
		d_end(0),
		evalue(DBL_MAX),
		bit_score(0.0),
		corrected_bit_score(0.0),
		approx_id(0.0),
#ifdef WITH_DNA
        mapping_quality(0),
        n_anchors(0),
#endif
		matrix(nullptr)
	{}

#ifdef WITH_DNA
    Hsp(const IntermediateRecord &r, unsigned query_source_len, Loc qlen, Loc tlen, const OutputFormat* output_format, const Stats::Blastn_Score *dna_score_builder = nullptr);
#else
	Hsp(const IntermediateRecord& r, unsigned query_source_len, Loc qlen, Loc tlen, const OutputFormat* output_format);
#endif
	Hsp(const ApproxHsp& h, Loc qlen, Loc tlen);

	struct Iterator
	{
		Iterator(const Hsp &parent) :
			query_pos(parent.query_range.begin_, Frame(parent.frame)),
			subject_pos(parent.subject_range.begin_),
			ptr_(parent.transcript.ptr()),
			count_(ptr_->count())
		{ }
		bool good() const
		{
			return *ptr_ != PackedOperation::terminator();
		}
		Iterator& operator++()
		{
			switch (op()) {
			case op_deletion:
				++subject_pos;
				break;
			case op_insertion:
				++query_pos;
				break;
			case op_match:
			case op_substitution:
				++query_pos;
				++subject_pos;
				break;
			case op_frameshift_forward:
				query_pos.shift_forward();
				break;
			case op_frameshift_reverse:
				query_pos.shift_back();
				break;
			}
			--count_;
			if (count_ == 0) {
				++ptr_;
				count_ = ptr_->count();
			}
			return *this;
		}
		EditOperation op() const
		{
			return ptr_->op();
		}
		bool to_next_match()
		{
			do {
				this->operator++();
			} while (good() && (op() == op_deletion || op() == op_insertion));
			return good();
		}
		TranslatedPosition query_pos;
		unsigned subject_pos;
	protected:
		const PackedOperation *ptr_;
		unsigned count_;
	};

	Iterator begin() const
	{
		return Iterator(*this);
	}

	Interval oriented_range() const
	{
		if (frame < 3)
			return Interval(query_source_range.begin_, query_source_range.end_ - 1);
		else
			return Interval(query_source_range.end_ - 1, query_source_range.begin_);
	}

	void set_translated_query_begin(unsigned oriented_query_begin, unsigned dna_len)
	{
		int f = frame <= 2 ? frame + 1 : 2 - frame;
		if (f > 0)
			query_range.begin_ = (oriented_query_begin - (f - 1)) / 3;
		else
			query_range.begin_ = (dna_len + f - oriented_query_begin) / 3;
	}

	void set_translated_query_end(unsigned oriented_query_end, unsigned dna_len)
	{
		int f = frame <= 2 ? frame + 1 : 2 - frame;
		if (f > 0)
			query_range.end_ = (oriented_query_end - 2 - (f - 1)) / 3 + 1;
		else
			query_range.end_ = (dna_len + f - oriented_query_end) / 3 + 1;
	}

	int blast_query_frame() const
	{
		return align_mode.query_translated ? (frame <= 2 ? (int)frame + 1 : 2 - (int)frame) : 0;
	}

	bool operator<(const Hsp &rhs) const
	{
		return score > rhs.score || (score == rhs.score && (d_begin < rhs.d_begin || (d_begin == rhs.d_begin && query_source_range.begin_ < rhs.query_source_range.begin_)));
	}

	static bool cmp_evalue(const Hsp& a, const Hsp& b) {
		return a.evalue < b.evalue || (a.evalue == b.evalue && a.score > b.score);
	}

	double id_percent() const
	{
		return (double)identities * 100.0 / (double)length;
	}

	double query_cover_percent(unsigned query_source_len) const
	{
		return (double)query_source_range.length() * 100 / query_source_len;
	}

	double subject_cover_percent(unsigned subject_len) const
	{
		return (double)subject_range.length() * 100 / subject_len;
	}

	static bool cmp_query_pos(const Hsp &x, const Hsp &y)
	{
		return x.query_range.begin_ < y.query_range.begin_;
	}

	int partial_score(const Hsp &h) const
	{
		const double overlap = std::max(subject_range.overlap_factor(h.subject_range), query_source_range.overlap_factor(h.query_source_range));
		return int((1 - overlap)*score);
	}

	bool envelopes(const DiagonalSegmentT &d, int dna_len) const
	{
		return query_source_range.contains(d.query_absolute_range(dna_len)) || subject_range.contains(d.subject_range());
	}

	bool is_enveloped_by(const Hsp &hsp, double p) const;
	bool is_enveloped_by(std::list<Hsp>::const_iterator begin, std::list<Hsp>::const_iterator end, double p) const;
	bool query_range_enveloped_by(const Hsp& hsp, double p) const;
	bool query_range_enveloped_by(std::list<Hsp>::const_iterator begin, std::list<Hsp>::const_iterator end, double p) const;
	bool is_weakly_enveloped_by(std::list<Hsp>::const_iterator begin, std::list<Hsp>::const_iterator end, int cutoff) const;
	void push_back(const DiagonalSegmentT &d, const TranslatedSequence &query, const Sequence& subject, bool reversed);
	void push_match(Letter q, Letter s, bool positive);
	void push_gap(EditOperation op, int length, const Letter *subject);
	void splice(const DiagonalSegmentT &d0, const DiagonalSegmentT &d1, const TranslatedSequence &query, const Sequence& subject, bool reversed);
	void set_begin(const DiagonalSegmentT &d, int dna_len);
	void set_end(const DiagonalSegmentT &d, int dna_len);
	void set_begin(int i, int j, Frame frame, int dna_len);
	void set_end(int i, int j, Frame frame, int dna_len);
	void clear();
	double approx_id_percent(const Sequence& query, const Sequence& target) const;
	bool is_identity(const Sequence& query, const Sequence& target) const;
	std::pair<Loc, Loc> min_range_len(double qcov, double tcov, Loc qlen, Loc tlen) const;

	bool is_weakly_enveloped(const Hsp &j) const;
	std::pair<int, int> diagonal_bounds() const;
	bool backtraced;
	int score, frame, length, identities, mismatches, positives, gap_openings, gaps, swipe_target, swipe_bin, d_begin, d_end;
#ifdef DP_STAT
	int reserved1, reserved2;
#endif
#if WITH_DNA
    int mapping_quality, n_anchors;
#endif
	Interval query_source_range, query_range, subject_source_range, subject_range;
	double evalue, bit_score, corrected_bit_score, approx_id;
	Sequence target_seq;
	const Stats::TargetMatrix* matrix;
	PackedTranscript transcript;
};

struct HspContext
{
	HspContext()
	{}
	HspContext(
		const Hsp& hsp,
		BlockId query_id,
		OId query_oid,
		const TranslatedSequence &query,
		const char *query_title,
		OId subject_oid,
		unsigned subject_len,
		const char* subject_title,
		int hit_num,
		int hsp_num,
		const Sequence &subject_seq,
		int ungapped_score = 0,
		const double query_self_aln_score = 0.0,
		const double target_self_aln_score = 0.0) :
		query(query),
		query_title(query_title),
		target_title(subject_title),
		query_id(query_id),
		query_oid(query_oid),
		subject_oid(subject_oid),
		query_len(query.source().length()),
		subject_len(subject_len),		
		hit_num(hit_num),
		hsp_num(hsp_num),
		query_self_aln_score(query_self_aln_score),
		target_self_aln_score(target_self_aln_score),
		subject_seq(subject_seq),
		hsp_(hsp)		
	{}
	struct Iterator : public Hsp::Iterator
	{
		Iterator(const HspContext &parent) :
			Hsp::Iterator(parent.hsp_),
			parent_(parent)
		{ }
		Letter query() const
		{
			return parent_.query[query_pos];
		}
		Letter subject() const
		{
			switch (op()) {
			case op_substitution:
			case op_deletion:
				return ptr_->letter();
			default:
				return query();
			}
		}
		char query_char() const
		{
			switch (op()) {
			case op_deletion:
				return '-';
			case op_frameshift_forward:
				return '\\';
			case op_frameshift_reverse:
				return '/';
			default:
				return value_traits.alphabet[(long)query()];
			}
		}
		char subject_char() const
		{
			switch (op()) {
			case op_insertion:
			case op_frameshift_forward:
			case op_frameshift_reverse:
				return '-';
			default:
				return value_traits.alphabet[(long)subject()];
			}
		}
		char midline_char(int score) const
		{
			switch (op()) {
			case op_match:
				return value_traits.alphabet[(long)query()];
			case op_substitution:
				return score > 0 ? '+' : ' ';
			default:
				return ' ';
			}
		}
	private:
		const HspContext &parent_;
	};
	Iterator begin() const
	{
		return Iterator(*this);
	}
	PackedTranscript::ConstIterator begin_old() const
	{
		return hsp_.transcript.begin();
	}
	unsigned score() const
	{ return hsp_.score; }
	double evalue() const
	{
		return hsp_.evalue;
	}
	double bit_score() const {
		return hsp_.bit_score;
	}
	double corrected_bit_score() const {
		return hsp_.corrected_bit_score;
	}
	double approx_id() const {
		return hsp_.approx_id;
	}
	unsigned frame() const
	{ return hsp_.frame; }
	unsigned length() const
	{ return hsp_.length; }
	unsigned identities() const
	{ return hsp_.identities; }
	unsigned mismatches() const
	{ return hsp_.mismatches; }
	unsigned positives() const
	{ return hsp_.positives; }
	unsigned gap_openings() const
	{ return hsp_.gap_openings; }
	unsigned gaps() const
	{ return hsp_.gaps; }
	double approx_id_percent() const {
		return hsp_.approx_id_percent(query.index(hsp_.frame), subject_seq);
	}
#if WITH_DNA
    unsigned mapping_quality() const
    { return hsp_.mapping_quality; }
    unsigned n_anchors() const
    { return hsp_.n_anchors; }
#endif
	double qcovhsp() const;
	double scovhsp() const;
	double id_percent() const;
	const Interval& query_source_range() const
	{ return hsp_.query_source_range; }
    const Interval& subject_source_range() const
    { return hsp_.subject_source_range; }
	const Interval& query_range() const
	{ return hsp_.query_range; }
	const Interval& subject_range() const
	{ return hsp_.subject_range; }
	Interval oriented_query_range() const
	{ return hsp_.oriented_range(); }
	int blast_query_frame() const
	{ return hsp_.blast_query_frame(); }
	PackedTranscript transcript() const
	{ return hsp_.transcript; }
	bool operator<(const HspContext& h) const {
		return query_oid < h.query_oid;
	}
#ifdef DP_STAT
	int reserved1() const {
		return hsp_.reserved1;
	}
	int reserved2() const {
		return hsp_.reserved2;
	}
#endif
	Hsp hsp() const {
		return hsp_;
	}
	HspContext& parse(const OutputFormat* output_format);

	TranslatedSequence query;
	std::string query_title, target_title;
	BlockId query_id;
	OId query_oid, subject_oid;
	Loc query_len, subject_len;
	unsigned hit_num, hsp_num;
	double query_self_aln_score, target_self_aln_score;
	Sequence subject_seq;
private:
	Hsp hsp_;
	friend HspContext deserialize(InputFile*);
};

inline void serialize(const HspContext& h, TextBuffer& buf) {
	buf.write(h.query_id);
	buf.write(h.query_oid);
	buf.write(h.subject_oid);
	buf.write_c_str(h.query_title.c_str());
	buf.write_c_str(h.target_title.c_str());
	buf.write(h.query_len);
	buf.write(h.subject_len);
	buf.write(h.identities());
	buf.write(h.mismatches());
	buf.write(h.positives());
	buf.write(h.gaps());
	buf.write(h.length());
	buf.write(h.gap_openings());
	buf.write(h.query_range().begin_);
	buf.write(h.query_range().end_);
	buf.write(h.subject_range().begin_);
	buf.write(h.subject_range().end_);
	buf.write(h.bit_score());
	buf.write(h.evalue());
	buf.write(h.score());
	buf.write(h.approx_id());
#if WITH_DNA
	buf.write(h.mapping_quality());
	buf.write(h.n_anchors());
#endif
}

inline HspContext deserialize(InputFile* file_) {
	HspContext h;
	file_->read(&h.query_id);
	file_->read(&h.query_oid);
	file_->read(&h.subject_oid);
	*file_ >> h.query_title;
	*file_ >> h.target_title;
	file_->read(&h.query_len);
	file_->read(&h.subject_len);
	file_->read(&h.hsp_.identities);
	file_->read(&h.hsp_.mismatches);
	file_->read(&h.hsp_.positives);
	file_->read(&h.hsp_.gaps);
	file_->read(&h.hsp_.length);
	file_->read(&h.hsp_.gap_openings);
	file_->read(&h.hsp_.query_range.begin_);
	file_->read(&h.hsp_.query_range.end_);
	file_->read(&h.hsp_.subject_range.begin_);
	file_->read(&h.hsp_.subject_range.end_);
	file_->read(&h.hsp_.bit_score);
	file_->read(&h.hsp_.evalue);
	file_->read(&h.hsp_.score);
	file_->read(&h.hsp_.approx_id);
#ifdef WITH_DNA
	file_->read(&h.hsp_.mapping_quality);
	file_->read(&h.hsp_.n_anchors);
#endif
	h.hsp_.query_source_range = h.hsp_.query_range;
	h.hsp_.subject_source_range = h.hsp_.subject_range;
	return h;
}