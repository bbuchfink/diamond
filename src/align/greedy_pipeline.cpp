/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#include "align.h"
#include "query_mapper.h"
#include "../util/map.h"
#include "../data/reference.h"

namespace ExtensionPipeline { namespace Greedy {

const bool log_ga = false;

struct Target : public ::Target
{

	void ungapped_stage(QueryMapper &mapper)
	{
		if (config.log_subject)
			cout << "Subject = " << ref_ids::get()[subject_id].c_str() << endl;
		std::stable_sort(mapper.seed_hits.begin() + begin, mapper.seed_hits.begin() + end, Seed_hit::compare_diag);
		typedef Map<vector<Seed_hit>::const_iterator, Seed_hit::Frame> Hit_map;
		Hit_map hit_map(mapper.seed_hits.begin() + begin, mapper.seed_hits.begin() + end);
		for (Hit_map::Iterator it = hit_map.begin(); it.valid(); ++it) {
			const unsigned frame = it.begin()->frame_;
			int score = greedy_align(mapper.query_seq(frame), mapper.profile[frame], mapper.query_cb[frame], subject, it.begin(), it.end(), log_ga, hsps, ts, frame);
			//target.filter_score = std::max(target.filter_score, (unsigned)target.traits[frame].score);
			if (score > filter_score) {
				filter_score = score;
			}
		}
		//target.filter_time = time;
	}

	void greedy_stage(QueryMapper &mapper, Statistics &stat, int cutoff)
	{
		if (config.log_subject)
			cout << "Subject = " << ref_ids::get()[subject_id].c_str() << endl;
#ifdef ENABLE_TIMING
		High_res_timer timer;
#endif
		filter_score = 0;
		ts.sort(Hsp_traits::cmp_diag);
		typedef Map<list<Hsp_traits>::const_iterator, Hsp_traits::Frame> Hsp_map;
		Hsp_map hsp_traits(ts.begin(), ts.end());
		list<Hsp_traits> t_out;
		hsps.clear();
		for (Hsp_map::Iterator it = hsp_traits.begin(); it.valid(); ++it) {
			const unsigned frame = it.begin()->frame;
			filter_score = std::max(filter_score, greedy_align(mapper.query_seq(frame), mapper.profile[frame], mapper.query_cb[frame], subject, log_ga, hsps, it.begin(), it.end(), t_out, cutoff, frame));
		}
		ts.clear();
		ts.splice(ts.begin(), t_out);
		//stat.inc(Statistics::TIME_GREEDY_EXT, timer.nanoseconds());
#ifdef ENABLE_TIMING
		target.filter_time = (float)timer.get();
#endif
	}

	void align_target(QueryMapper &mapper, Statistics &stat)
	{
		const size_t n = end - begin;
		if (config.log_subject)
			cout << "Subject = " << ref_ids::get()[subject_id].c_str() << endl;

		stat.inc(Statistics::CELLS, mapper.query_seq(0).length() * subject.length());

		if (filter_score == 0)
			return;
		if (config.ext == Config::more_greedy)
			hsps.emplace_back(filter_score);
		else {
			const int qlen = (int)mapper.query_seq(0).length(),
				band_plus = qlen <= 50 ? 0 : 16;
			hsps.clear();
			for (list<Hsp_traits>::const_iterator i = ts.begin(); i != ts.end(); ++i) {
				if (log_ga) {
					cout << "i_begin=" << i->query_range.begin_ << " j_begin=" << i->subject_range.begin_ << " d_min=" << i->d_min << " d_max=" << i->d_max << endl;
				}
				hsps.emplace_back();
				hsps.back().frame = i->frame;
				banded_sw(mapper.query_seq(i->frame), subject, i->d_min - band_plus, i->d_max + band_plus + 1, 0, (int)subject.length(), hsps.back());
				if (config.comp_based_stats) {
					const int score = (int)hsps.back().score + mapper.query_cb[i->frame](hsps.back());
					hsps.back().score = (unsigned)std::max(0, score);
				}
			}
		}

		if (!hsps.empty())
			stat.inc(Statistics::OUT_HITS);

		for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end(); ++i)
			for (list<Hsp>::iterator j = hsps.begin(); j != hsps.end();)
				if (j != i && j->is_weakly_enveloped(*i)) {
					stat.inc(Statistics::ERASED_HITS);
					j = hsps.erase(j);
				}
				else
					++j;

		//const float time = (float)timer.getElapsedTimeInMicroSec() + target.filter_time;

		for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end(); ++i) {
			i->time = filter_time;
			i->query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(i->query_range.begin_, Frame(i->frame)), TranslatedPosition(i->query_range.end_, Frame(i->frame)), mapper.source_query_len);
		}

		hsps.sort();
		if (hsps.size() > 0)
			filter_score = hsps.front().score;

		ts.clear();
		for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end(); ++i)
			ts.emplace_back(i->query_source_range);

		if (config.use_smith_waterman && !hsps.empty()) {
			int score;
			for (unsigned f = 0; f < align_mode.query_contexts; ++f) {
				needleman_wunsch(mapper.query_seq(f), subject, score, Local(), int());
				hsps.front().sw_score = std::max((unsigned)score, hsps.front().sw_score);
			}
			stat.inc(Statistics::SQUARED_ERROR, (stat_type)pow((int)hsps.front().sw_score - (int)hsps.front().score, 2));
		}
	}

};

Target& Pipeline::target(size_t i)
{
	return (Target&)(this->targets[i]);
}

void Pipeline::run(Statistics &stat)
{
	if (n_targets() == 0)
		return;
	stat.inc(Statistics::TARGET_HITS0, n_targets());
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).ungapped_stage(*this);
	if (config.ext == Config::most_greedy)
		return;
	fill_source_ranges();
	rank_targets(config.rank_ratio == -1 ? (query_seq(0).length() > 50 ? 0.6 : 0.9) : config.rank_ratio, config.rank_factor == -1 ? 1e3 : config.rank_factor);
	stat.inc(Statistics::TARGET_HITS1, n_targets());
	const int cutoff = int(raw_score_cutoff() * config.score_ratio);
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).greedy_stage(*this, stat, cutoff);
	fill_source_ranges();
	rank_targets(config.rank_ratio2 == -1 ? (query_seq(0).length() > 50 ? 0.95 : 1.0) : config.rank_ratio2, config.rank_factor == -1 ? 1e3 : config.rank_factor);
	stat.inc(Statistics::TARGET_HITS2, n_targets());
	for (size_t i = 0; i < n_targets(); ++i)
		target(i).align_target(*this, stat);
}

}}