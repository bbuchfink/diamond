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

#ifndef TARGET_CULLING_H_
#define TARGET_CULLING_H_

#include <vector>
#include "../align/query_mapper.h"

struct Target_culling
{
	Target_culling():
		top_score_(0)
	{}
	Target_culling(const Target *first):
		n_(0),
		top_score_(first ? first->filter_score : 0)
	{}
	virtual int cull(const Target &t)
	{
		return config.output_range((unsigned)n_, t.filter_score, top_score_) ? use : finished;
	}
	virtual void add(Target *t)
	{
		++n_;
	}
	enum { finished = 0, next = 1, use = 2};
	static Target_culling* get(const Target* first);
private:
	size_t n_;
	const int top_score_;
};

struct Overlap_culling : Target_culling
{
	Overlap_culling()
	{}
	virtual int cull(const Target &t)
	{
		return t.is_enveloped(targets_.begin(), targets_.end(), config.query_overlap_culling / 100.0, int(double(t.filter_score) / (1.0 - config.toppercent / 100.0)))
			? next : use;
	}
	virtual void add(Target *t)
	{
		targets_.push_back(t);
	}
	virtual void add(unsigned n, const Hsp_data &hsp)
	{
		if (n >= targets_.size()) {
			targets_.push_back(new Target(hsp.score));
		}
		targets_[n]->ts.push_back(Hsp_traits(hsp));
	}
private:
	std::vector<Target*> targets_;
};

#endif