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
		std::stable_sort(hits, hits_end);
		filter_score = hits[0].ungapped.score;
	}

	void add(QueryMapper &mapper, vector<DpTarget> &vf, vector<DpTarget> &vr)
	{
		vector<Seed_hit>::iterator hits = mapper.seed_hits.begin() + begin, hits_end = mapper.seed_hits.begin() + end;
		const int d = hits[0].diagonal();
		const DpTarget t(subject, std::max(d - 32, -int(subject.length() - 1)), std::min(d + 32, int(mapper.query_seq(0).length() - 1)), &hsps, subject_id);
		if (hits[0].frame_ < 3)
			vf.push_back(t);
		else
			vr.push_back(t);
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

void Pipeline::run(Statistics &stat)
{
	//cout << "Query=" << query_ids::get()[this->query_id].c_str() << endl;
	vector<DpTarget> vf, vr;
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).ungapped_stage(*this);
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).add(*this, vf, vr);
	std::sort(vf.begin(), vf.end());
	std::sort(vr.begin(), vr.end());
	for (vector<DpTarget>::iterator i = vf.begin(); i < vf.end(); i += 8) {
		banded_3frame_swipe(translated_query, FORWARD, i, std::min(i + 8, vf.end()));
	}
	for (vector<DpTarget>::iterator i = vr.begin(); i < vr.end(); i += 8)
		banded_3frame_swipe(translated_query, REVERSE, i, std::min(i + 8, vr.end()));
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).finish(*this);
}

}}
