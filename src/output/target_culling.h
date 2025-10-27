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
#include <vector>
#include <set>
#include <map>
#include "align/legacy/query_mapper.h"
#include "util/geo/interval_partition.h"
#include "output.h"

struct TargetCulling
{
	virtual std::pair<int, double> cull(const Target &t) const = 0;
	virtual int cull(const std::vector<IntermediateRecord> &target_hsp, const std::set<TaxId> &taxon_ids) const = 0;
	virtual void add(const Target &t) = 0;
	virtual void add(const std::vector<IntermediateRecord> &target_hsp, const std::set<TaxId> &taxon_ids) = 0;
	virtual ~TargetCulling() = default;
	enum { FINISHED = 0, NEXT = 1, INCLUDE = 2};
	static TargetCulling* get(const int64_t max_target_seqs);
};

struct GlobalCulling : public TargetCulling
{
	GlobalCulling(const int64_t max_target_seqs) :
		max_target_seqs_(max_target_seqs),
		n_(0),
		top_score_(0)
	{}
	virtual std::pair<int, double> cull(const Target &t) const
	{
		if (top_score_ == 0)
			return { INCLUDE, 0 };
		if (config.taxon_k) {
			unsigned taxons_exceeded = 0;
			for (unsigned i : t.taxon_rank_ids) {
				auto it = taxon_count_.find(i);
				if (it != taxon_count_.end() && it->second >= config.taxon_k)
					++taxons_exceeded;
			}
			if (taxons_exceeded == t.taxon_rank_ids.size())
				return { NEXT,0 };
		}
		if (config.toppercent.present())
			return { (1.0 - score_matrix.bitscore(t.filter_score) / top_score_) * 100.0 <= config.toppercent ? INCLUDE : FINISHED, 0 };
		else
			return { n_ < max_target_seqs_ ? INCLUDE : FINISHED,0 };
	}
	virtual int cull(const std::vector<IntermediateRecord> &target_hsp, const std::set<TaxId> &taxon_ids) const
	{
		if (top_score_ == 0.0)
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
		if (config.global_ranking_targets)
			return n_ < config.global_ranking_targets ? INCLUDE : FINISHED;
		else if (config.toppercent.present())
			return (1.0 - score_matrix.bitscore(target_hsp[0].score) / top_score_) * 100.0 <= config.toppercent ? INCLUDE : FINISHED;
		else
			return n_ < max_target_seqs_ ? INCLUDE : FINISHED;
	}
	virtual void add(const Target &t)
	{
		if (top_score_ == 0)
			top_score_ = score_matrix.bitscore(t.filter_score);
		++n_;
		if (config.taxon_k)
			for (unsigned i : t.taxon_rank_ids)
				++taxon_count_[i];
	}
	virtual void add(const std::vector<IntermediateRecord> &target_hsp, const std::set<TaxId> &taxon_ids)
	{
		if (top_score_ == 0)
			top_score_ = score_matrix.bitscore(target_hsp[0].score);
		++n_;
		if (config.taxon_k)
			for (unsigned i : taxon_ids)
				++taxon_count_[i];
	}
	virtual ~GlobalCulling() = default;
private:
	const int64_t max_target_seqs_;
	int64_t n_;
	double top_score_;
	std::map<unsigned, unsigned> taxon_count_;
};

struct RangeCulling : public TargetCulling
{
	RangeCulling(const int64_t max_target_seqs) :
		p_(max_target_seqs)
	{}
	virtual std::pair<int, double> cull(const Target &t) const
	{
		int c = 0, l = 0;
		for (std::list<Hsp>::const_iterator i = t.hsps.begin(); i != t.hsps.end(); ++i) {
			if (config.toppercent.blank()) {
				c += p_.covered(i->query_source_range);
			}
			else {
				const int cutoff = int((double)i->score / (1.0 - config.toppercent / 100.0));
				c += p_.covered(i->query_source_range, cutoff, IntervalPartition::MaxScore());
			}
			l += i->query_source_range.length();
		}
		const double cov = (double)c / l;
		return std::make_pair(cov * 100.0 < config.query_range_cover ? INCLUDE : NEXT, cov);
	}
	virtual int cull(const std::vector<IntermediateRecord> &target_hsp, const std::set<TaxId> &taxon_ids) const
	{
		int c = 0, l = 0;
		for (std::vector<IntermediateRecord>::const_iterator i = target_hsp.begin(); i != target_hsp.end(); ++i) {
			if (config.toppercent.blank())
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
	virtual void add(const std::vector<IntermediateRecord> &target_hsp, const std::set<TaxId> &taxon_ids)
	{
		for (std::vector<IntermediateRecord>::const_iterator i = target_hsp.begin(); i != target_hsp.end(); ++i)
			p_.insert(i->absolute_query_range(), i->score);
	}
	virtual ~RangeCulling() = default;
private:
	IntervalPartition p_;
};