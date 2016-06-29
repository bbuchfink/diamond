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

#ifndef ALIGN_H_
#define ALIGN_H_

#include <vector>
#include "../search/trace_pt_buffer.h"
#include "../util/task_queue.h"
#include "../basic/statistics.h"

using std::vector;

struct local_match : public Hsp_data
{
	typedef vector<local_match>::iterator iterator;
	local_match() :
		query_anchor_(0),
		subject_(0)
	{ }
	local_match(int score) :
		Hsp_data(score)
	{ }
	local_match(int query_anchor, int subject_anchor, const Letter *subject, unsigned total_subject_len = 0) :
		total_subject_len_(total_subject_len),
		query_anchor_(query_anchor),
		subject_anchor(subject_anchor),
		subject_(subject)
	{ }
	local_match(unsigned len, unsigned query_begin, unsigned query_len, unsigned subject_len, unsigned gap_openings, unsigned identities, unsigned mismatches, signed subject_begin, signed score) :
		Hsp_data(score),
		query_anchor_(0),
		subject_(0)
	{ }
	void merge(const local_match &right, const local_match &left);
	bool pass_through(const Diagonal_segment &d);
	bool is_weakly_enveloped(const local_match &j);
	unsigned total_subject_len_;
	signed query_anchor_, subject_anchor;
	const Letter *subject_;
};

struct Segment
{
	Segment(int score,
		unsigned frame,
		local_match *traceback = 0,
		unsigned subject_id = std::numeric_limits<unsigned>::max()) :
		score_(score),
		frame_(frame),
		traceback_(traceback),
		subject_id_(subject_id),
		next_(0),
		top_score_(0)
	{ }
	Strand strand() const
	{
		return frame_ < 3 ? FORWARD : REVERSE;
	}
	bool operator<(const Segment &rhs) const
	{
		return top_score_ > rhs.top_score_
			|| (top_score_ == rhs.top_score_
				&& (subject_id_ < rhs.subject_id_ || (subject_id_ == rhs.subject_id_ && (score_ > rhs.score_ || (score_ == rhs.score_ && traceback_->score > rhs.traceback_->score)))));
	}
	static bool comp_subject(const Segment& lhs, const Segment &rhs)
	{
		return lhs.subject_id_ < rhs.subject_id_ || (lhs.subject_id_ == rhs.subject_id_ && lhs.score_ > rhs.score_);
	}
	struct Subject
	{
		unsigned operator()(const Segment& x) const
		{
			return x.subject_id_;
		}
	};
	int						score_;
	unsigned				frame_;
	local_match			   *traceback_;
	unsigned				subject_id_;
	Segment				   *next_;
	int						top_score_;
};

void align_sequence_simple(vector<Segment> &matches,
	Statistics &stat,
	vector<local_match> &local,
	unsigned *padding,
	size_t db_letters,
	unsigned dna_len,
	Trace_pt_buffer::Vector::iterator &begin,
	Trace_pt_buffer::Vector::iterator &end);

void align_sequence_anchored(vector<Segment> &matches,
	Statistics &stat,
	vector<local_match> &local,
	unsigned *padding,
	size_t db_letters,
	unsigned dna_len,
	Trace_pt_buffer::Vector::iterator &begin,
	Trace_pt_buffer::Vector::iterator &end);

struct Output_writer
{
	Output_writer(Output_stream* f) :
		f_(f)
	{ }
	void operator()(Text_buffer &buf)
	{
		f_->write(buf.get_begin(), buf.size());
		buf.clear();
	}
private:
	Output_stream* const f_;
};

template<typename _buffer>
struct Ring_buffer_sink
{
	Ring_buffer_sink(Output_stream *output_file):
		writer(output_file),
		queue(config.threads_ * 32, writer)
	{}
	bool get(size_t &i, _buffer *& buffer, Trace_pt_list::Query_range &query_range)
	{
		return queue.get(i, buffer, query_range);
	}
	void push(size_t i)
	{
		queue.push(i);
	}
private:
	Output_writer writer;
	Task_queue<_buffer, Output_writer> queue;
};

#endif
