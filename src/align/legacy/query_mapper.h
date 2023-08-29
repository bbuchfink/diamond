/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include <queue>
#include <vector>
#include <list>
#include <set>
#include <float.h>
#include "../../data/queries.h"
#include "../../util/ptr_vector.h"
#include "../../dp/dp.h"
#include "../../data/reference.h"
#include "../../util/hsp/approx_hsp.h"
#include "../../basic/match.h"
#include "../../run/config.h"
#include "../search/hit.h"

struct Seed_hit
{
	Seed_hit()
	{}
	Seed_hit(unsigned frame, unsigned subject, unsigned subject_pos, unsigned query_pos, const DiagonalSegment &ungapped) :
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
		const DiagonalSegmentT d(ungapped, ::Frame(frame_));
		for (std::list<Hsp>::const_iterator i = begin; i != end; ++i)
			if (i->envelopes(d, dna_len))
				return true;
		return false;
	}
	DiagonalSegmentT diagonal_segment() const
	{
		return DiagonalSegmentT(ungapped, ::Frame(frame_));
	}
	Interval query_source_range(int dna_len) const
	{
		return diagonal_segment().query_absolute_range(dna_len);
	}
	Strand strand() const
	{
		return ::Frame(frame_).strand;
	}
	static bool compare_pos(const Seed_hit &x, const Seed_hit &y)
	{
		return DiagonalSegment::cmp_subject_end(x.ungapped, y.ungapped);
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
	DiagonalSegment ungapped;
	unsigned prefix_score;
};

struct Target
{
	Target(int filter_score, double filter_evalue):
		filter_score(filter_score),
		filter_evalue(filter_evalue)
	{}
	Target(size_t begin, unsigned subject_id, const Sequence& subject, const std::set<TaxId> &taxon_rank_ids) :
		subject_block_id(subject_id),
		subject(subject),
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
		for (std::list<ApproxHsp>::iterator i = ts.begin(); i != ts.end(); ++i)
			i->query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(i->query_range.begin_, Frame(i->frame)), TranslatedPosition(i->query_range.end_, Frame(i->frame)), (int)query_len);
	}
	void add_ranges(std::vector<int32_t> &v);
	bool is_outranked(const std::vector<int32_t> &v, double treshold);
	bool envelopes(const ApproxHsp &t, double p) const;
	bool is_enveloped(const Target &t, double p) const;
	bool is_enveloped(PtrVector<Target>::const_iterator begin, PtrVector<Target>::const_iterator end, double p, int min_score) const;
	void inner_culling();
	void apply_filters(int dna_len, int subject_len, const char *query_title);
	unsigned subject_block_id;
	Sequence subject;
	int filter_score;
	double filter_evalue;
	float filter_time;
	bool outranked;
	size_t begin, end;
	std::list<Hsp> hsps;
	std::list<ApproxHsp> ts;
	Seed_hit top_hit;
	std::set<TaxId> taxon_rank_ids;

	enum { INTERVAL = 64 };
};

struct QueryMapper
{
	QueryMapper(size_t query_id, Search::Hit* begin, Search::Hit* end, const Search::Config &metadata);
	void init();
	bool generate_output(TextBuffer &buffer, Statistics &stat, const Search::Config& cfg);
	void rank_targets(double ratio, double factor, const int64_t max_target_seqs);
	void score_only_culling(const int64_t max_target_seqs);
	int64_t n_targets() const
	{
		return (int64_t)targets.size();
	}
	bool finished() const
	{
		return targets_finished == targets.size();
	}
	Sequence query_seq(unsigned frame) const
	{
		return metadata.query->seqs()[(size_t)query_id*(size_t)align_mode.query_contexts + frame];
	}
	void fill_source_ranges()
	{
		for (size_t i = 0; i < targets.size(); ++i)
			targets[i].fill_source_ranges(source_query_len);
	}
	virtual void run(Statistics &stat, const Search::Config& cfg) = 0;
	virtual ~QueryMapper() {}

	std::pair<Search::Hit*, Search::Hit*> source_hits;
	unsigned query_id, targets_finished, next_target;
	unsigned source_query_len, unaligned_from;
	PtrVector<Target> targets;
	std::vector<Seed_hit> seed_hits;
	std::vector<Bias_correction> query_cb;
	TranslatedSequence translated_query;
	const Search::Config &metadata;
	bool target_parallel;

private:

	static std::pair<Search::Hit*, Search::Hit*> get_query_data();
	unsigned count_targets();
	Sequence query_source_seq() const
	{
		return align_mode.query_translated ? metadata.query->source_seqs()[query_id] : metadata.query->seqs()[query_id];
	}
	void load_targets();
	
};
