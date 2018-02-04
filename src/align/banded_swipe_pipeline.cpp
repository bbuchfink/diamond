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

	interval ungapped_query_range(int query_dna_len) const
	{
		const int i0 = std::max((int)top_hit.query_pos_ - (int)top_hit.subject_pos_, 0),
			i1 = std::min((int)top_hit.query_pos_ + (int)subject.length() - (int)top_hit.subject_pos_, query_dna_len / 3);
		const Frame f = Frame(top_hit.frame_);
		return TranslatedPosition::absolute_interval(TranslatedPosition(i0, f), TranslatedPosition(i1, f), query_dna_len);
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
	}

	void reset()
	{
		hsps.clear();
	}

	void finish(QueryMapper &mapper)
	{
		inner_culling(mapper.raw_score_cutoff());
	}

	bool is_outranked(const IntervalPartition &ip, int source_query_len, double rr) const
	{
		const interval r = ungapped_query_range(source_query_len);
		if (config.toppercent == 100.0) {
			const int min_score = int((double)filter_score / rr);
			return (double)ip.covered(r, min_score, IntervalPartition::MinScore()) / r.length() * 100.0 >= config.query_range_cover;
		}
		else {
			const int min_score = int((double)filter_score / rr / (1.0 - config.toppercent / 100.0));
			return (double)ip.covered(r, min_score, IntervalPartition::MaxScore()) / r.length() * 100.0 >= config.query_range_cover;
		}
	}

};

void Pipeline::range_ranking()
{
	const double rr = config.rank_ratio == -1 ? 0.4 : config.rank_ratio;
	std::stable_sort(targets.begin(), targets.end(), Target::compare);
	IntervalPartition ip((int)std::min(config.max_alignments, (uint64_t)INT_MAX));
	for (PtrVector< ::Target>::iterator i = targets.begin(); i < targets.end();) {
		Target* t = ((Target*)*i);
		if (t->is_outranked(ip, source_query_len, rr)) {
			if (config.benchmark_ranking) {
				t->outranked = true;
				++i;
			}
			else
				i = targets.erase(i, i + 1);
		}
		else {
			ip.insert(t->ungapped_query_range(source_query_len), t->filter_score);
			++i;
		}			
	}
}

Target& Pipeline::target(size_t i)
{
	return (Target&)(this->targets[i]);
}

void Pipeline::run_swipe(bool score_only)
{
	vector<DpTarget> vf, vr;
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).add(*this, vf, vr);
	banded_3frame_swipe(translated_query, FORWARD, vf.begin(), vf.end(), this->dp_stat, score_only);
	banded_3frame_swipe(translated_query, REVERSE, vr.begin(), vr.end(), this->dp_stat, score_only);
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
	if (!config.query_range_culling)
		rank_targets(config.rank_ratio == -1 ? 0.4 : config.rank_ratio, config.rank_factor == -1.0 ? 1e3 : config.rank_factor);
	else
		range_ranking();

	if (n_targets() > config.max_alignments) {
		stat.inc(Statistics::TARGET_HITS1, n_targets());
		run_swipe(true);
		for (size_t i = 0; i < n_targets(); ++i)
			target(i).set_filter_score();
		score_only_culling();
	}

	stat.inc(Statistics::TARGET_HITS2, n_targets());
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).reset();
	run_swipe(false);
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).finish(*this);
}

}}
