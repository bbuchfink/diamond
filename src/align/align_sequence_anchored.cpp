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

// #define ENABLE_LOGGING_AS

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
		subject_pos_(subject_pos),
		query_pos_(query_pos),
		extensions_ (0),
		ungapped(ungapped_extension(subject, subject_pos, query_pos, query)),
		hsp_(hsp)
	{ }
	friend std::ostream& operator<<(std::ostream &os, const local_trace_point &p)
	{
		cout << "(subject_pos=" << p.subject_pos_ << ",query_pos_=" << p.query_pos_ << ")";
		return os;
	}
	int diagonal() const
	{
		return (int)subject_pos_ - (int)query_pos_;
	}
	bool operator<(const local_trace_point &rhs) const
	{
		return ungapped.score > rhs.ungapped.score;
	}
	unsigned subject_pos_, query_pos_, extensions_;
	Diagonal_segment ungapped;
	local_match *hsp_;
};

struct Subject_seq
{
	Subject_seq(unsigned id, unsigned begin) :
		id(id),
		begin(begin),
		filter_score (0)
	{}
	bool operator<(const Subject_seq &rhs) const
	{
		return filter_score > rhs.filter_score;
	}
	unsigned id, begin, end, filter_score;
};

void load_local_trace_points(vector<local_trace_point> &v, vector<Subject_seq> &subjects, Trace_pt_buffer::Vector::iterator &begin, Trace_pt_buffer::Vector::iterator &end, const sequence &query)
{
	size_t subject_id = std::numeric_limits<size_t>::max();
	for (Trace_pt_buffer::Vector::iterator i = begin; i != end; ++i) {
		std::pair<size_t, size_t> l = ref_seqs::data_->local_position(i->subject_);
		if (l.first != subject_id) {
			if (v.size() > 0)
				subjects.back().end = (unsigned)v.size();
			subject_id = l.first;
			subjects.push_back(Subject_seq((unsigned)subject_id, (unsigned)v.size()));
		}
		v.push_back(local_trace_point((unsigned)l.first, (unsigned)l.second, i->seed_offset_, query));
		subjects.back().filter_score += v.back().ungapped.score;
	}
	subjects.back().end = (unsigned)v.size();
}

void rank_subjects(vector<Subject_seq> &subjects, vector<local_trace_point> &tp)
{
	std::sort(subjects.begin(), subjects.end());
	int score = (int)((double)subjects[std::min(subjects.size(), config.max_alignments) - 1].filter_score * config.rank_ratio);
	unsigned i = 0;
	for (; i < subjects.size(); ++i)
		if (subjects[i].filter_score < score)
			break;

	subjects.erase(subjects.begin() + std::min((unsigned)subjects.size(), std::max((unsigned)(config.max_alignments*config.rank_factor), i)), subjects.end());

	for (vector<Subject_seq>::iterator i = subjects.begin(); i != subjects.end(); ++i)
		std::sort(tp.begin() + i->begin, tp.begin() + i->end);
}

bool include(vector<local_trace_point>::iterator begin,
	vector<local_trace_point>::iterator end,
	const local_trace_point &p,
	unsigned band)
{
	if (config.extend_all)
		return true;
	for (vector<local_trace_point>::iterator i = begin; i<end; ++i) {
		if (i->hsp_ != 0
			//&& (i->hsp_->score_ == 0 || i->hsp_->pass_through(Diagonal_segment(p.query_pos_, p.subject_pos_, 1), transcript_buf)))
			&& (i->hsp_->score == 0
				|| p.ungapped.is_enveloped(i->ungapped)
				|| i->hsp_->pass_through(p.ungapped)))
			return false;
	}
	return true;
}

void load_subject_seqs(unsigned subject_id,
	vector<local_match> &dst,
	vector<local_trace_point>::iterator begin,
	vector<local_trace_point>::iterator end,
	const sequence &query,
	unsigned query_len,
	unsigned band,
	Statistics &stat)
{
	//if (begin->extensions_ >= sqrt(query_len)/2)
	/*if(begin->extensions_ >= ref_seqs::get().length(begin->subject_) * 10)
		return;*/
	double sim = std::max((double)(end - begin) / query_len * ref_seqs::data_->length(subject_id),
		(double)(end - begin) * query_len / ref_seqs::data_->length(subject_id));

	if (sim > 200) {
		if (begin->hsp_)
			return;
		//cout << end - begin << ' ' << query_len << ' ' << sim << endl;
		stat.inc(Statistics::HIGH_SIM);
		dst.push_back(local_match(begin->query_pos_, begin->subject_pos_, ref_seqs::data_->data(ref_seqs::data_->position(subject_id, begin->subject_pos_))));
		begin->hsp_ = &dst.back();
		return;
	}
	/*if (end - begin > 100) {
		cout << "n=" << end - begin << " querylen=" << query_len << " reflen=" << ref_seqs::data_->length(subject_id) << " sim=" << (double)(end - begin) / query_len * ref_seqs::data_->length(subject_id) << endl;
	}*/
	for (vector<local_trace_point>::iterator i = begin; i < end; ++i) {
		if (i->hsp_ == 0 && include(begin, end, *i, band)) {
			dst.push_back(local_match(i->query_pos_, i->subject_pos_, ref_seqs::data_->data(ref_seqs::data_->position(subject_id, i->subject_pos_))));
#ifdef ENABLE_LOGGING_AS
			cout << i->query_pos_ << ' ' << i->subject_ << ' ' << i->subject_pos_ << endl;
			cout << ref_ids::get()[i->subject_].c_str() << endl;
#endif
			i->hsp_ = &dst.back();
			//++begin->extensions_;
			//begin->extensions_ += i->hsp_->len_;
		}
	}
}

void load_subject_seqs(const vector<Subject_seq> &subjects,
	vector<local_match> &dst,
	vector<local_trace_point> &src,
	const sequence &query,
	unsigned query_len,
	unsigned band,
	Statistics &stat)
{
	for (vector<Subject_seq>::const_iterator i = subjects.begin(); i != subjects.end();++i)
		load_subject_seqs(i->id, dst, src.begin() + i->begin, src.begin() + i->end, query, query_len, band, stat);
}

void align_sequence_anchored(vector<Segment> &matches,
	Statistics &stat,
	vector<local_match> &local,
	unsigned *padding,
	size_t db_letters,
	unsigned dna_len,
	Trace_pt_buffer::Vector::iterator &begin,
	Trace_pt_buffer::Vector::iterator &end)
{
	static TLS_PTR vector<local_trace_point> *trace_points;
	static TLS_PTR vector<Subject_seq> *subject_ptr;

	vector<local_trace_point> &trace_pt(get_tls(trace_points));
	vector<Subject_seq> &subjects(get_tls(subject_ptr));

	trace_pt.clear();
	subjects.clear();

	std::sort(begin, end, hit::cmp_subject);
	
	const unsigned q_num(begin->query_);
	const sequence query(query_seqs::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = (unsigned)query.length();
	padding[frame] = config.read_padding(query_len);
	unsigned aligned = 0;
	uint64_t cell_updates = 0;
	load_local_trace_points(trace_pt, subjects, begin, end, query);
	rank_subjects(subjects, trace_pt);

	while (true) {
		size_t local_begin = local.size();
		load_subject_seqs(subjects, local, trace_pt, query, query_len, padding[frame], stat);
		
		if (local.size() - local_begin == 0) {
			stat.inc(Statistics::OUT_HITS, aligned);
			stat.inc(Statistics::DUPLICATES, trace_pt.size() - aligned);
			break;
		}
		aligned += (unsigned)(local.size() - local_begin);

		for (size_t i = local_begin; i < local.size(); ++i) {
			floating_sw(&query[local[i].query_anchor_],
				local[i],
				padding[frame],
				score_matrix.rawscore(config.gapped_xdrop),
				config.gap_open + config.gap_extend,
				config.gap_extend,
				cell_updates,
				Traceback());
		}
	}

	for (vector<Subject_seq>::const_iterator s = subjects.begin(); s != subjects.end();++s)
		for (vector<local_trace_point>::iterator i = trace_pt.begin() + s->begin; i != trace_pt.begin() + s->end; ++i)
			if (i->hsp_) {
				for (vector<local_trace_point>::iterator j = trace_pt.begin() + s->begin; j != trace_pt.begin() + s->end; ++j)
					if (i != j && j->hsp_ && i->hsp_->is_weakly_enveloped(*j->hsp_)) {
						i->hsp_ = 0;
						break;
					}
			}

	for (vector<Subject_seq>::const_iterator s = subjects.begin(); s != subjects.end(); ++s)
		for (vector<local_trace_point>::iterator i = trace_pt.begin() + s->begin; i != trace_pt.begin() + s->end; ++i)
			if (i->hsp_) {
				i->hsp_->set_source_range(frame, dna_len);
				matches.push_back(Segment(i->hsp_->score, frame, i->hsp_, s->id));
			}
}