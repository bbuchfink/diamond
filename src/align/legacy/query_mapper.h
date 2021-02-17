/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef QUERY_MAPPER_H_
#define QUERY_MAPPER_H_

#include <queue>
#include <vector>
#include <list>
#include <set>
#include <float.h>
#include "../../search/trace_pt_buffer.h"
#include "../../data/queries.h"
#include "../../util/ptr_vector.h"
#include "../../dp/dp.h"
#include "../../data/reference.h"
#include "../../basic/parameters.h"
#include "../../data/metadata.h"
#include "../../dp/hsp_traits.h"
#include "../../basic/match.h"

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
	bool is_enveloped(std::list<Hsp>::const_iterator begin, std::list<Hsp>::const_iterator end, int dna_len) const
	{
		const DiagonalSegment d(ungapped, ::Frame(frame_));
		for (std::list<Hsp>::const_iterator i = begin; i != end; ++i)
			if (i->envelopes(d, dna_len))
				return true;
		return false;
	}
	DiagonalSegment diagonal_segment() const
	{
		return DiagonalSegment(ungapped, ::Frame(frame_));
	}
	interval query_source_range(int dna_len) const
	{
		return diagonal_segment().query_absolute_range(dna_len);
	}
	Strand strand() const
	{
		return ::Frame(frame_).strand;
	}
	static bool compare_pos(const Seed_hit &x, const Seed_hit &y)
	{
		return Diagonal_segment::cmp_subject_end(x.ungapped, y.ungapped);
	}
	static bool compare_diag(const Seed_hit &x, const Seed_hit &y)
	{
		return x.frame_ < y.frame_ || (x.frame_ == y.frame_ && (x.diagonal() < y.diagonal() || (x.diagonal() == y.diagonal() && x.ungapped.j < y.ungapped.j)));
	}
	static bool compare_diag_strand(const Seed_hit &x, const Seed_hit &y)
	{
		return x.strand() < y.strand() || (x.strand() == y.strand() && (x.diagonal() < y.diagonal() || (x.diagonal() == y.diagonal() && x.ungapped.j < y.ungapped.j)));
	}
	static bool compare_diag_strand2(const Seed_hit &x, const Seed_hit &y)
	{
		return x.strand() < y.strand() || (x.strand() == y.strand() && (x.diagonal() < y.diagonal() || (x.diagonal() == y.diagonal() && x.subject_pos_ < y.subject_pos_)));
	}
	struct Frame
	{
		unsigned operator()(const Seed_hit &x) const
		{
			return x.frame_;
		}
	};

	unsigned frame_, subject_, subject_pos_, query_pos_;
	Diagonal_segment ungapped;
	unsigned prefix_score;
};

struct Target
{
	Target(int filter_score, double filter_evalue):
		filter_score(filter_score),
		filter_evalue(filter_evalue)
	{}
	Target(size_t begin, unsigned subject_id, const std::set<unsigned> &taxon_rank_ids) :
		subject_block_id(subject_id),
		subject(ref_seqs::get()[subject_id]),
		filter_score(0),
		filter_evalue(DBL_MAX),
		outranked(false),
		begin(begin),
		taxon_rank_ids(taxon_rank_ids)
	{}
	static bool compare_evalue(Target* lhs, Target *rhs)
	{
		return lhs->filter_evalue < rhs->filter_evalue || (lhs->filter_evalue == rhs->filter_evalue && compare_score(lhs, rhs));
	}
	static bool compare_score(Target* lhs, Target* rhs)
	{
		return lhs->filter_score > rhs->filter_score || (lhs->filter_score == rhs->filter_score && lhs->subject_block_id < rhs->subject_block_id);
	}
	void fill_source_ranges(size_t query_len)
	{
		for (list<Hsp_traits>::iterator i = ts.begin(); i != ts.end(); ++i)
			i->query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(i->query_range.begin_, Frame(i->frame)), TranslatedPosition(i->query_range.end_, Frame(i->frame)), (int)query_len);
	}
	void add_ranges(vector<int32_t> &v);
	bool is_outranked(const vector<int32_t> &v, double treshold);
	bool envelopes(const Hsp_traits &t, double p) const;
	bool is_enveloped(const Target &t, double p) const;
	bool is_enveloped(PtrVector<Target>::const_iterator begin, PtrVector<Target>::const_iterator end, double p, int min_score) const;
	void inner_culling();
	void apply_filters(int dna_len, int subject_len, const char *query_title, const char *ref_title);
	unsigned subject_block_id;
	Sequence subject;
	int filter_score;
	double filter_evalue;
	float filter_time;
	bool outranked;
	size_t begin, end;
	list<Hsp> hsps;
	list<Hsp_traits> ts;
	Seed_hit top_hit;
	std::set<unsigned> taxon_rank_ids;

	enum { INTERVAL = 64 };
};

struct QueryMapper
{
	QueryMapper(const Parameters &params, size_t query_id, hit* begin, hit* end, const Metadata &metadata, bool target_parallel = false);
	void init();
	bool generate_output(TextBuffer &buffer, Statistics &stat);
	void rank_targets(double ratio, double factor);
	void score_only_culling();
	size_t n_targets() const
	{
		return targets.size();
	}
	bool finished() const
	{
		return targets_finished == targets.size();
	}
	Sequence query_seq(unsigned frame) const
	{
		return query_seqs::get()[query_id*align_mode.query_contexts + frame];
	}
	void fill_source_ranges()
	{
		for (size_t i = 0; i < targets.size(); ++i)
			targets[i].fill_source_ranges(source_query_len);
	}
	virtual void run(Statistics &stat) = 0;
	virtual ~QueryMapper() {}

	const Parameters &parameters;
	pair<hit*, hit*> source_hits;
	unsigned query_id, targets_finished, next_target;
	unsigned source_query_len, unaligned_from;
	PtrVector<Target> targets;
	vector<Seed_hit> seed_hits;
	vector<Bias_correction> query_cb;
	TranslatedSequence translated_query;
	bool target_parallel;
	const Metadata &metadata;

private:

	static pair<hit*, hit*> get_query_data();
	unsigned count_targets();
	Sequence query_source_seq() const
	{
		return align_mode.query_translated ? query_source_seqs::get()[query_id] : query_seqs::get()[query_id];
	}
	void load_targets();
	
};

#endif