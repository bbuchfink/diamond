/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
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
#include <vector>
#include <set>
#include <map>
#include "../align/legacy/query_mapper.h"
#include "../util/interval_partition.h"
#include "output.h"

struct TargetCulling
{
	virtual int cull(const Target &t) const = 0;
	virtual int cull(const vector<IntermediateRecord> &target_hsp, const std::set<unsigned> &taxon_ids) const = 0;
	virtual void add(const Target &t) = 0;
	virtual void add(const vector<IntermediateRecord> &target_hsp, const std::set<unsigned> &taxon_ids) = 0;
	virtual ~TargetCulling() = default;
	enum { FINISHED = 0, NEXT = 1, INCLUDE = 2};
	static TargetCulling* get();
};

struct GlobalCulling : public TargetCulling
{
	GlobalCulling() :
		n_(0),
		top_score_(0)
	{}
	virtual int cull(const Target &t) const
	{
		if (top_score_ == 0)
			return INCLUDE;
		if (config.taxon_k) {
			unsigned taxons_exceeded = 0;
			for (unsigned i : t.taxon_rank_ids) {
				auto it = taxon_count_.find(i);
				if (it != taxon_count_.end() && it->second >= config.taxon_k)
					++taxons_exceeded;
			}
			if (taxons_exceeded == t.taxon_rank_ids.size())
				return NEXT;
		}
		if (config.toppercent < 100.0)
			return (1.0 - (double)t.filter_score / top_score_) * 100.0 <= config.toppercent ? INCLUDE : FINISHED;
		else
			return n_ < config.max_alignments ? INCLUDE : FINISHED;
	}
	virtual int cull(const vector<IntermediateRecord> &target_hsp, const std::set<unsigned> &taxon_ids) const
	{
		if (top_score_ == 0)
			return INCLUDE;
		if (config.taxon_k) {
			unsigned taxons_exceeded = 0;
			for (unsigned i : taxon_ids) {
				auto it = taxon_count_.find(i);
				if (it != taxon_count_.end() && it->second >= config.taxon_k)
					++taxons_exceeded;
			}
			if (taxons_exceeded == taxon_ids.size())
				return NEXT;
		}
		if (config.toppercent < 100.0)
			return (1.0 - (double)target_hsp[0].score / top_score_) * 100.0 <= config.toppercent ? INCLUDE : FINISHED;
		else
			return n_ < config.max_alignments ? INCLUDE : FINISHED;
	}
	virtual void add(const Target &t)
	{
		if (top_score_ == 0)
			top_score_ = t.filter_score;
		++n_;
		if (config.taxon_k)
			for (unsigned i : t.taxon_rank_ids)
				++taxon_count_[i];
	}
	virtual void add(const vector<IntermediateRecord> &target_hsp, const std::set<unsigned> &taxon_ids)
	{
		if (top_score_ == 0)
			top_score_ = target_hsp[0].score;
		++n_;
		if (config.taxon_k)
			for (unsigned i : taxon_ids)
				++taxon_count_[i];
	}
	virtual ~GlobalCulling() = default;
private:
	size_t n_;
	int top_score_;
	std::map<unsigned, unsigned> taxon_count_;
};

struct RangeCulling : public TargetCulling
{
	RangeCulling() :
		p_((int)std::min(config.max_alignments, (uint64_t)INT_MAX))
	{}
	virtual int cull(const Target &t) const
	{
		int c = 0, l = 0;
		for (std::list<Hsp>::const_iterator i = t.hsps.begin(); i != t.hsps.end(); ++i) {
			if (config.toppercent == 100.0) {
				c += p_.covered(i->query_source_range);
			}
			else {
				const int cutoff = int((double)i->score / (1.0 - config.toppercent / 100.0));
				c += p_.covered(i->query_source_range, cutoff, IntervalPartition::MaxScore());
			}
			l += i->query_source_range.length();
		}
		return (double)c / l * 100.0 < config.query_range_cover ? INCLUDE : NEXT;
	}
	virtual int cull(const std::vector<IntermediateRecord> &target_hsp, const std::set<unsigned> &taxon_ids) const
	{
		int c = 0, l = 0;
		for (std::vector<IntermediateRecord>::const_iterator i = target_hsp.begin(); i != target_hsp.end(); ++i) {
			if (config.toppercent == 100.0)
				c += p_.covered(i->absolute_query_range());
			else {
				const int cutoff = int((double)i->score / (1.0 - config.toppercent / 100.0));
				c += p_.covered(i->absolute_query_range(), cutoff, IntervalPartition::MaxScore());
			}
			l += i->absolute_query_range().length();
		}
		return (double)c / l * 100.0 < config.query_range_cover ? INCLUDE : NEXT;
	}
	virtual void add(const Target &t)
	{
		for (std::list<Hsp>::const_iterator i = t.hsps.begin(); i != t.hsps.end(); ++i)
			p_.insert(i->query_source_range, i->score);
	}
	virtual void add(const vector<IntermediateRecord> &target_hsp, const std::set<unsigned> &taxon_ids)
	{
		for (std::vector<IntermediateRecord>::const_iterator i = target_hsp.begin(); i != target_hsp.end(); ++i)
			p_.insert(i->absolute_query_range(), i->score);
	}
	virtual ~RangeCulling() = default;
private:
	IntervalPartition p_;
};