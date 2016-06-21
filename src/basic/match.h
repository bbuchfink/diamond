/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
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

#ifndef MATCH_H_
#define MATCH_H_

#include <limits>
#include "sequence.h"
#include "../util/async_buffer.h"
#include "edit_transcript.h"
#include "packed_loc.h"
#include "../util/system.h"
#include "value.h"

enum Strand { FORWARD, REVERSE };

inline interval normalized_range(unsigned pos, int len, Strand strand)
{
	return strand == FORWARD
			? interval (pos, pos + len)
			: interval (pos + 1 + len, pos + 1);
}

struct Diagonal_segment
{
	Diagonal_segment(int query_pos, int subject_pos, int len, int score):
		query_pos(query_pos),
		subject_pos(subject_pos),
		len(len),
		score (score)
	{}
	interval query_range() const
	{
		return interval(query_pos, query_pos + len);
	}
	interval subject_range() const
	{
		return interval(subject_pos, subject_pos + len);
	}
	int diag() const
	{
		return subject_pos - query_pos;
	}
	bool is_enveloped(const Diagonal_segment &x) const
	{
		return score <= x.score
			&& query_range().overlap_factor(x.query_range()) == 1
			&& subject_range().overlap_factor(x.subject_range()) == 1;
	}
	int query_pos, subject_pos, len, score;
};

#pragma pack(1)

struct hit
{
	typedef uint32_t Seed_offset;

	unsigned	query_;
	Packed_loc	subject_;
	Seed_offset	seed_offset_;
	hit() :
		query_(),
		subject_(),
		seed_offset_()
	{ }
	hit(unsigned query, Packed_loc subject, Seed_offset seed_offset) :
		query_(query),
		subject_(subject),
		seed_offset_(seed_offset)
	{ }
	bool operator<(const hit &rhs) const
	{
		return query_ < rhs.query_;
	}
	bool blank() const
	{
		return subject_ == 0;
	}
	unsigned operator%(unsigned i) const
	{
		return (query_ / query_contexts()) % i;
	}
	unsigned operator/(size_t i) const
	{
		return (query_ / query_contexts()) / (unsigned)i;
	}
	int64_t global_diagonal() const
	{
		return (int64_t)subject_ - (int64_t)seed_offset_;
	}
	template<unsigned _d>
	static unsigned query_id(const hit& x)
	{
		return x.query_ / _d;
	}
	template<unsigned _d>
	struct Query_id
	{
		unsigned operator()(const hit& x) const
		{
			return query_id<_d>(x);
		}
	};
	static bool cmp_subject(const hit &lhs, const hit &rhs)
	{
		return lhs.subject_ < rhs.subject_
			|| (lhs.subject_ == rhs.subject_ && lhs.seed_offset_ < rhs.seed_offset_);
	}
	static bool cmp_normalized_subject(const hit &lhs, const hit &rhs)
	{
		const uint64_t x = (uint64_t)lhs.subject_ + (uint64_t)rhs.seed_offset_, y = (uint64_t)rhs.subject_ + (uint64_t)lhs.seed_offset_;
		return x < y || (x == y && lhs.seed_offset_ < rhs.seed_offset_);
	}
	friend std::ostream& operator<<(std::ostream &s, const hit &me)
	{
		s << me.query_ << '\t' << me.subject_ << '\t' << me.seed_offset_ << '\n';
		return s;
	}
} PACKED_ATTRIBUTE ;

#pragma pack()

struct Hsp_data
{
	unsigned score, frame, length, identities, mismatches, positives, gap_openings, gaps;
	interval query_source_range, query_range, subject_range;
	Packed_transcript transcript;
};

struct local_match
{
	typedef vector<local_match>::iterator iterator;
	local_match():
		len_ (),
		query_begin_ (),
		subject_len_ (),
		gap_openings_ (),
		identities_ (),
		mismatches_ (),
		subject_begin_ (),
		score_ (),
		query_len_ (),
		query_anchor_ (0),
		subject_ (0)
	{ }
	local_match(int score):
		len_ (0),
		query_begin_ (0),
		subject_len_ (0),
		gap_openings_ (0),
		identities_ (0),
		mismatches_ (0),
		subject_begin_ (0),
		score_ (score),
		query_len_ (0),
		query_anchor_ (0),
		subject_ (0)
	{ }
	local_match(int query_anchor, int subject_anchor, const Letter *subject, unsigned total_subject_len = 0):
		len_ (0),
		query_begin_ (0),
		subject_len_ (0),
		gap_openings_ (0),
		identities_ (0),
		mismatches_ (0),
		total_subject_len_(total_subject_len),
		subject_begin_ (0),
		score_ (0),
		query_len_ (0),
		query_anchor_ (query_anchor),
		subject_anchor (subject_anchor),
		subject_ (subject)
	{ }
	local_match(unsigned len, unsigned query_begin, unsigned query_len, unsigned subject_len, unsigned gap_openings, unsigned identities, unsigned mismatches, signed subject_begin, signed score):
		len_ (len),
		query_begin_ (query_begin),
		subject_len_ (subject_len),
		gap_openings_ (gap_openings),
		identities_ (identities),
		mismatches_ (mismatches),
		subject_begin_ (subject_begin),
		score_ (score),
		query_len_ (query_len),
		query_anchor_ (0),
		subject_ (0)
	{ }
	local_match& operator+=(const local_match& rhs)
	{
		add(rhs);
		transcript_right_ = rhs.transcript_right_;
		return *this;
	}
	local_match& operator-=(const local_match& rhs)
	{
		add(rhs);
		query_begin_ = rhs.query_len_;
		subject_begin_ = rhs.subject_len_;
		transcript_left_ = rhs.transcript_right_;
		return *this;
	}
	void add(const local_match &rhs)
	{
		len_ += rhs.len_;
		subject_len_ += rhs.subject_len_;
		gap_openings_ += rhs.gap_openings_;
		identities_ += rhs.identities_;
		mismatches_ += rhs.mismatches_;
		score_ += rhs.score_;
		query_len_ += rhs.query_len_;
	}
	interval query_range(Strand strand) const
	{ return normalized_range(query_begin_, query_len_, strand); }
	interval subject_range() const
	{ return normalized_range(subject_begin_, subject_len_, FORWARD); }
	void print(const sequence &query, const sequence &subject, const vector<char> &buf) const
	{
		std::cout << "Score = " << score_ << std::endl;
		::print(std::cout, &query[query_begin_], &subject[subject_begin_], transcript_right_, transcript_left_, buf);
	}
	bool pass_through(const Diagonal_segment &d, const vector<char> &transcript_buf);
	bool is_weakly_enveloped(const local_match &j);
	unsigned len_, query_begin_, subject_len_, gap_openings_, identities_, mismatches_, total_subject_len_;
	signed subject_begin_, score_, query_len_, query_anchor_, subject_anchor;
	const Letter *subject_;
	Edit_transcript transcript_right_, transcript_left_;
};

struct Segment
{
	Segment(int score,
			unsigned frame,
			local_match *traceback = 0,
			unsigned subject_id = std::numeric_limits<unsigned>::max()):
		score_ (score),
		frame_ (frame),
		traceback_ (traceback),
		subject_id_ (subject_id),
		next_ (0),
		top_score_ (0)
	{ }
	Strand strand() const
	{ return frame_ < 3 ? FORWARD : REVERSE; }
	interval query_range() const
	{ return traceback_->query_range(strand()); }
	interval subject_range() const
	{ return traceback_->subject_range(); }
	bool operator<(const Segment &rhs) const
	{ return top_score_ > rhs.top_score_
			|| (top_score_ == rhs.top_score_
			&& (subject_id_ < rhs.subject_id_ || (subject_id_ == rhs.subject_id_ && (score_ > rhs.score_ || (score_ == rhs.score_ && traceback_->score_ > rhs.traceback_->score_))))); }
	static bool comp_subject(const Segment& lhs, const Segment &rhs)
	{ return lhs.subject_id_ < rhs.subject_id_ || (lhs.subject_id_ == rhs.subject_id_ && lhs.score_ > rhs.score_); }
	struct Subject
	{
		unsigned operator()(const Segment& x) const
		{ return x.subject_id_; }
	};
	int						score_;
	unsigned				frame_;
	local_match			   *traceback_;
	unsigned				subject_id_;
	Segment				   *next_;
	int						top_score_;
};

#endif /* MATCH_H_ */
