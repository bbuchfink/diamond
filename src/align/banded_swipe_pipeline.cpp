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

#include "align.h"
#include "../dp/dp.h"
#include "../util/interval_partition.h"

namespace ExtensionPipeline { namespace BandedSwipe {

struct Target : public ::Target
{

	void ungapped_stage(QueryMapper &mapper)
	{
		vector<Seed_hit>::iterator hits = mapper.seed_hits.begin() + begin, hits_end = mapper.seed_hits.begin() + end;
		top_hit = hits[0];
		for (vector<Seed_hit>::iterator i = hits + 1; i < hits_end; ++i)
			if (i->ungapped.score > top_hit.ungapped.score)
				top_hit = *i;
		filter_score = top_hit.ungapped.score;
	}

	void add_strand(QueryMapper &mapper, vector<DpTarget> &v, vector<Seed_hit>::iterator begin, vector<Seed_hit>::iterator end)
	{
		if (end == begin) return;
		const int band = config.padding, d_min = -int(subject.length() - 1), d_max = int(mapper.query_seq(0).length() - 1);
		vector<Seed_hit>::iterator i = begin;
		int d0 = std::max(i->diagonal() - band, d_min), d1 = std::min(i->diagonal() + band, d_max);
		++i;
		for (; i < end; ++i) {
			if (i->diagonal() - d1 <= band) {
				d1 = std::min(i->diagonal() + band, d_max);
			}
			else {
				v.push_back(DpTarget(subject, d0, d1, &hsps, subject_id));
				d0 = std::max(i->diagonal() - band, d_min);
				d1 = std::min(i->diagonal() + band, d_max);
			}
		}
		v.push_back(DpTarget(subject, d0, d1, &hsps, subject_id));
	}

	void add(QueryMapper &mapper, vector<DpTarget> &vf, vector<DpTarget> &vr)
	{
		vector<Seed_hit>::iterator hits = mapper.seed_hits.begin() + begin, hits_end = mapper.seed_hits.begin() + end;
		Strand strand = top_hit.strand();
		std::stable_sort(hits, hits_end, Seed_hit::compare_diag_strand);
		for (vector<Seed_hit>::iterator i = hits; i < hits_end; ++i)
			if (i->strand() == REVERSE) {
				if (strand == FORWARD)
					add_strand(mapper, vf, hits, i);
				else
					add_strand(mapper, vr, i, hits_end);
				return;
			}
		if (strand == FORWARD) add_strand(mapper, vf, hits, hits_end);
		//const int d = hits[0].diagonal();
		//const DpTarget t(subject, std::max(d - 32, -int(subject.length() - 1)), std::min(d + 32, int(mapper.query_seq(0).length() - 1)), &hsps, subject_id);		
	}

	void set_filter_score()
	{
		filter_score = 0;
		for (list<Hsp>::const_iterator i = hsps.begin(); i != hsps.end(); ++i)
			filter_score = std::max(filter_score, (int)i->score);
		hsps.clear();
	}

	void finish(QueryMapper &mapper)
	{
		inner_culling(mapper.raw_score_cutoff());
	}

};

Target& Pipeline::target(size_t i)
{
	return (Target&)(this->targets[i]);
}

void Pipeline::run_swipe(bool score_only)
{
	vector<DpTarget> vf, vr;
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).add(*this, vf, vr);
	std::sort(vf.begin(), vf.end());
	std::sort(vr.begin(), vr.end());
	if (!vf.empty())
		for (vector<DpTarget>::iterator i = vf.begin(); i < vf.end(); i += 8) {
			banded_3frame_swipe(translated_query, FORWARD, i, std::min(i + 8, vf.end()), this->dp_stat, score_only);
		}
	if (!vr.empty())
		for (vector<DpTarget>::iterator i = vr.begin(); i < vr.end(); i += 8)
			banded_3frame_swipe(translated_query, REVERSE, i, std::min(i + 8, vr.end()), this->dp_stat, score_only);
}

void Pipeline::run(Statistics &stat)
{
	//cout << "Query=" << query_ids::get()[this->query_id].c_str() << endl;
	Config::set_option(config.padding, 32);
	if (n_targets() == 0)
		return;
	stat.inc(Statistics::TARGET_HITS0, n_targets());
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).ungapped_stage(*this);
	run_swipe(true);
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).set_filter_score();
	score_only_culling();
	stat.inc(Statistics::TARGET_HITS1, n_targets());
	run_swipe(false);
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).finish(*this);
}

}}
