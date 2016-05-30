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

#include <stddef.h>
#include "../dp/floating_sw.h"
#include "../data/queries.h"
#include "align.h"
#include "../data/reference.h"
#include "match_func.h"
#include "../util/map.h"
#include "../dp/dp.h"

using std::vector;

Diagonal_segment ungapped_extension(unsigned subject, unsigned subject_pos, unsigned query_pos, const sequence &query)
{
	const Letter* s = ref_seqs::data_->data(ref_seqs::data_->position(subject, subject_pos)),
		*q = &query[query_pos];
	unsigned delta, len;
	int score = xdrop_ungapped(q, s, delta, len);
	return Diagonal_segment(query_pos - delta, subject_pos - delta, len, score);
}

struct local_trace_point
{
	local_trace_point(unsigned subject, unsigned subject_pos, unsigned query_pos, const sequence& query, local_match *hsp = 0) :
		subject_(subject),
		subject_pos_(subject_pos),
		query_pos_(query_pos),
		ungapped(ungapped_extension(subject, subject_pos, query_pos, query)),
		hsp_(hsp)
	{ }
	struct Subject
	{
		unsigned operator()(const local_trace_point &x) const
		{
			return x.subject_;
		}
	};
	friend std::ostream& operator<<(std::ostream &os, const local_trace_point &p)
	{
		cout << "(subject=" << p.subject_ << ",subject_pos=" << p.subject_pos_ << ",query_pos_=" << p.query_pos_ << ")";
		return os;
	}
	int diagonal() const
	{
		return (int)subject_pos_ - (int)query_pos_;
	}
	unsigned subject_, subject_pos_, query_pos_;
	Diagonal_segment ungapped;
	local_match *hsp_;
};

void load_local_trace_points(vector<local_trace_point> &v, Trace_pt_buffer::Vector::iterator &begin, Trace_pt_buffer::Vector::iterator &end, const sequence &query)
{
	for (Trace_pt_buffer::Vector::iterator i = begin; i != end; ++i) {
		std::pair<size_t, size_t> l = ref_seqs::data_->local_position(i->subject_);
		v.push_back(local_trace_point((unsigned)l.first, (unsigned)l.second, i->seed_offset_, query));
		//cout << v.back() << endl;
	}
}

bool include(vector<local_trace_point>::iterator begin,
	vector<local_trace_point>::iterator end,
	const local_trace_point &p,
	unsigned band,
	const vector<char> &transcript_buf)
{
	if (config.extend_all)
		return true;
	for (vector<local_trace_point>::iterator i = begin; i<end; ++i) {
		if (i->hsp_ != 0
			//&& (i->hsp_->score_ == 0 || i->hsp_->pass_through(Diagonal_segment(p.query_pos_, p.subject_pos_, 1), transcript_buf)))
			&& (i->hsp_->score_ == 0
				|| p.ungapped.is_enveloped(i->ungapped)
				|| i->hsp_->pass_through(p.ungapped, transcript_buf)))
			return false;
	}
	return true;
}

void load_subject_seqs(vector<local_match> &dst,
	vector<local_trace_point>::iterator begin,
	vector<local_trace_point>::iterator end,
	const sequence &query,
	unsigned query_len,
	unsigned band,
	const vector<char> &transcript_buf)
{
	for (vector<local_trace_point>::iterator i = begin; i < end; ++i) {
		if (i->hsp_ == 0 && include(begin, end, *i, band, transcript_buf)) {
			dst.push_back(local_match(i->query_pos_, i->subject_pos_, ref_seqs::data_->data(ref_seqs::data_->position(i->subject_, i->subject_pos_))));
			i->hsp_ = &dst.back();
			//cout << *i << " i\n";
		}
		else {
			//cout << *i << " x\n";
		}
	}
}

void load_subject_seqs(vector<local_match> &dst,
	vector<local_trace_point> &src,
	const sequence &query,
	unsigned query_len,
	unsigned band,
	const vector<char> &transcript_buf)
{
	typedef Map<vector<local_trace_point>::iterator, local_trace_point::Subject> trace_pt_map;
	trace_pt_map trace_pt(src.begin(), src.end());
	trace_pt_map::Iterator it = trace_pt.begin();
	while (it.valid()) {
		load_subject_seqs(dst, it.begin(), it.end(), query, query_len, band, transcript_buf);
		++it;
	}
}
void align_sequence_anchored(vector<Segment> &matches,
	Statistics &stat,
	vector<local_match> &local,
	unsigned *padding,
	size_t db_letters,
	unsigned dna_len,
	Trace_pt_buffer::Vector::iterator &begin,
	Trace_pt_buffer::Vector::iterator &end,
	vector<char> &transcript_buf)
{
	static TLS_PTR vector<local_trace_point> *trace_points;

	vector<local_trace_point> &trace_pt(get_tls(trace_points));

	std::sort(begin, end, hit::cmp_subject);
	trace_pt.clear();
	const unsigned q_num(begin->query_);
	const sequence query(query_seqs::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = (unsigned)query.length();
	padding[frame] = config.read_padding(query_len);
	unsigned aligned = 0;
	load_local_trace_points(trace_pt, begin, end, query);

	while (true) {
		const local_match::iterator local_begin = local.end();
		load_subject_seqs(local, trace_pt, query, query_len, padding[frame], transcript_buf);
		
		if (local.end() - local_begin == 0) {
			stat.inc(Statistics::OUT_HITS, aligned);
			stat.inc(Statistics::DUPLICATES, trace_pt.size() - aligned);
			break;
		}
		aligned += (unsigned)(local.end() - local_begin);

		for (vector<local_match>::iterator i = local_begin; i < local.end(); ++i) {
			floating_sw(&query[i->query_anchor_],
				*i,
				padding[frame],
				score_matrix.rawscore(config.gapped_xdrop),
				config.gap_open + config.gap_extend,
				config.gap_extend,
				transcript_buf,
				Traceback());
			anchored_transform(*i, i->subject_anchor, i->query_anchor_);
		}
	}

	typedef Map<vector<local_trace_point>::iterator, local_trace_point::Subject> trace_pt_map;
	trace_pt_map map(trace_pt.begin(), trace_pt.end());
	trace_pt_map::Iterator it = map.begin();
	while (it.valid()) {
		for (vector<local_trace_point>::iterator i = it.begin(); i != it.end(); ++i)
			if (i->hsp_) {
				for (vector<local_trace_point>::iterator j = it.begin(); j != it.end(); ++j)
					if (i != j && j->hsp_ && i->hsp_->is_weakly_enveloped(*j->hsp_)) {
						i->hsp_ = 0;
						break;
					}
			}
		++it;
	}

	for (vector<local_trace_point>::iterator i = trace_pt.begin(); i != trace_pt.end(); ++i)
		if (i->hsp_) {
			to_source_space(*(i->hsp_), frame, dna_len);
			matches.push_back(Segment(i->hsp_->score_, frame, i->hsp_, i->subject_));
		}
}