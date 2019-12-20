/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

#ifndef TARGET_H_
#define TARGET_H_

#include <array>
#include <set>
#include <vector>
#include <stdint.h>
#include <list>
#include "../search/trace_pt_buffer.h"
#include "../basic/diagonal_segment.h"
#include "../basic/const.h"
#include "../dp/hsp_traits.h"
#include "../dp/comp_based_stats.h"

namespace Extension {

struct WorkTarget {
	WorkTarget(const sequence &seq) :
		seq(seq),
		filter_score(0),
		outranked(false)
	{}
	bool operator<(const WorkTarget &t) const {
		return filter_score > t.filter_score;
	}
	sequence seq;
	int filter_score;
	bool outranked;
	std::array<std::list<Hsp_traits>, MAX_CONTEXT> hsp;
};

std::vector<WorkTarget> ungapped_stage(const sequence *query_seq, const Bias_correction *query_cb, Trace_pt_list::iterator begin, Trace_pt_list::iterator end);
void rank_targets(std::vector<WorkTarget> &targets, double ratio, double factor);

struct Context {
};

struct Target {

	Target(const sequence &seq):
		filter_score(0),
		seq(seq)
	{}
	Target(Trace_pt_list::const_iterator begin, Trace_pt_list::const_iterator end, uint64_t target_offset, const sequence *query_seq, const sequence &target_seq, const std::set<unsigned> &taxon_rank_ids);
	void add_hit(std::list<Hsp> &list, std::list<Hsp>::iterator it) {
		std::list<Hsp> &l = hsp[it->frame];
		l.splice(l.end(), list, it);
		filter_score = std::max(filter_score, (int)l.back().score);
	}
	bool operator<(const Target &t) const {
		return filter_score > t.filter_score;
	}
	void inner_culling();

	sequence seq;
	int filter_score;
	std::array<std::list<Hsp>, MAX_CONTEXT> hsp;
	//const std::set<unsigned> taxon_rank_ids;

};

void score_only_culling(std::vector<Target> &targets);
std::vector<Target> align(const std::vector<WorkTarget> &targets, const sequence *query_seq, const Bias_correction *query_cb);
std::vector<Target> align(const std::vector<Target> &targets, const sequence *query_seq, const Bias_correction *query_cb);

}

#endif