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

struct Target
{
	Target(int score):
		filter_score(score)
	{}
	Target(size_t begin, unsigned subject_id) :
		subject_id(subject_id),
		filter_score(0),
		outranked(false),
		begin(begin)
	{}		
	static bool compare(Target* lhs, Target *rhs)
	{
		return lhs->filter_score > rhs->filter_score || (lhs->filter_score == rhs->filter_score && lhs->subject_id < rhs->subject_id);
	}
	void fill_source_ranges(size_t query_len)
	{
		for (list<Hsp_traits>::iterator i = ts.begin(); i != ts.end(); ++i)
			i->query_source_range = align_mode.query_translated ? untranslate_range(i->query_range, i->frame, query_len) : i->query_range;
	}
	bool envelopes(const Hsp_traits &t, double p) const;
	bool is_enveloped(const Target &t, double p) const;
	bool is_enveloped(Ptr_vector<Target>::const_iterator begin, Ptr_vector<Target>::const_iterator end, double p, int min_score) const;
	unsigned subject_id;
	int filter_score;
	float filter_time;
	bool outranked;
	size_t begin, end;
	list<Hsp_data> hsps;
	list<Hsp_traits> ts;
};

struct Query_mapper
{
	Query_mapper();
	Query_mapper(size_t query_id, Trace_pt_list::iterator begin, Trace_pt_list::iterator end);
	void init();
	void get_prefilter_score(size_t idx);
	void align_target(size_t idx, Statistics &stat);
	void align_targets(Statistics &stat);
	bool generate_output(Text_buffer &buffer, Statistics &stat);
	void ungapped_stage(size_t idx);
	void greedy_stage(size_t idx, Statistics &stat, int cutoff);
	void rank_targets(double ratio);
	int raw_score_cutoff() const;
	size_t n_targets() const
	{
		return targets.size();
	}
	bool finished() const
	{
		return targets_finished == targets.size();
	}
	sequence query_seq(unsigned frame) const
	{
		return query_seqs::get()[query_id*align_mode.query_contexts + frame];
	}
	void fill_source_ranges()
	{
		for (size_t i = 0; i < targets.size(); ++i)
			targets[i].fill_source_ranges(source_query_len);
	}
	pair<Trace_pt_list::iterator, Trace_pt_list::iterator> source_hits;
	unsigned query_id, targets_finished, next_target;
private:
	static pair<Trace_pt_list::iterator, Trace_pt_list::iterator> get_query_data();
	unsigned count_targets();
	sequence query_source_seq() const
	{
		return align_mode.query_translated ? query_source_seqs::get()[query_id] : query_seqs::get()[query_id];
	}
	void load_targets();

	unsigned source_query_len, unaligned_from;
	vector<Seed_hit> seed_hits;
	Ptr_vector<Target> targets;
	vector<Bias_correction> query_cb;
	vector<Long_score_profile> profile;
	
};

#endif