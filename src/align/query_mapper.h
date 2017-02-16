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

#ifndef QUERY_MAPPER_H_
#define QUERY_MAPPER_H_

#include <queue>
#include <vector>
#include <list>
#include "../util/tinythread.h"
#include "../search/trace_pt_buffer.h"
#include "../data/queries.h"
#include "../util/ptr_vector.h"
#include "../dp/dp.h"

using std::vector;
using std::pair;
using std::list;

struct Seed_hit
{
	Seed_hit()
	{}
	Seed_hit(unsigned frame, unsigned subject, unsigned subject_pos, unsigned query_pos, const Diagonal_segment &ungapped) :
		frame_(frame),
		subject_(subject),
		subject_pos_(subject_pos),
		query_pos_(query_pos),
		ungapped(ungapped),
		prefix_score(ungapped.score)
	{ }
	int diagonal() const
	{
		return (int)query_pos_ - (int)subject_pos_;
	}
	bool operator<(const Seed_hit &rhs) const
	{
		return ungapped.score > rhs.ungapped.score;
	}
	static bool compare_pos(const Seed_hit &x, const Seed_hit &y)
	{
		return Diagonal_segment::cmp_subject_end(x.ungapped, y.ungapped);
	}
	unsigned frame_, subject_, subject_pos_, query_pos_;
	Diagonal_segment ungapped;
	unsigned prefix_score;
};

struct Target
{
	Target(size_t begin, unsigned subject_id) :
		subject_id(subject_id),
		filter_score(0),
		begin(begin)
	{}		
	static bool compare(Target* lhs, Target *rhs)
	{
		return lhs->filter_score > rhs->filter_score;
	}
	unsigned subject_id, filter_score;
	size_t begin, end;
	list<Hsp_data> hsps;
};

struct Query_mapper
{
	Query_mapper();
	void init();
	void get_prefilter_score(size_t idx);
	void align_target(size_t idx, Statistics &stat);
	bool generate_output(Text_buffer &buffer, Statistics &stat);
	size_t n_targets() const
	{
		return targets.size();
	}
	bool finished() const
	{
		return targets_finished == targets.size();
	}
	pair<Trace_pt_list::iterator, Trace_pt_list::iterator> source_hits;
	unsigned query_id, targets_finished, next_target;
private:
	static pair<Trace_pt_list::iterator, Trace_pt_list::iterator> get_query_data();
	unsigned count_targets();
	sequence query_seq(unsigned frame) const
	{
		return query_seqs::get()[query_id*align_mode.query_contexts + frame];
	}
	sequence query_source_seq() const
	{
		return align_mode.query_translated ? query_source_seqs::get()[query_id] : query_seqs::get()[query_id];
	}
	void load_targets();
	void rank_targets();

	unsigned source_query_len, unaligned_from;
	vector<Seed_hit> seed_hits;
	Ptr_vector<Target> targets;
	vector<Bias_correction> query_cb;
	vector<Long_score_profile> profile;
	
};

#endif