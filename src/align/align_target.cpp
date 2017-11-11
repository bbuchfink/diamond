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

#include "query_mapper.h"
#include "../data/reference.h"
#include "../dp/floating_sw.h"
#include "../util/map.h"
#include "../util/high_res_timer.h"

using std::list;

// #define ENABLE_TIMING

const bool log_ga = false;

void Query_mapper::ungapped_stage(size_t idx)
{
	Target& target = targets[idx];
	const string subject_id(ref_ids::get()[target.subject_id].c_str());
	if (config.log_subject)
		cout << "Subject = " << subject_id << endl;
	std::stable_sort(seed_hits.begin() + target.begin, seed_hits.begin() + target.end, Seed_hit::compare_diag);
	typedef Map<vector<Seed_hit>::const_iterator, Seed_hit::Frame> Hit_map;
	Hit_map hit_map(seed_hits.begin() + target.begin, seed_hits.begin() + target.end);
	const sequence subject = ref_seqs::get()[target.subject_id];
	for (Hit_map::Iterator it = hit_map.begin(); it.valid(); ++it) {
		const unsigned frame = it.begin()->frame_;
		int score = greedy_align(query_seq(frame), profile[frame], query_cb[frame], subject, it.begin(), it.end(), log_ga, target.hsps, target.ts, frame);
		//target.filter_score = std::max(target.filter_score, (unsigned)target.traits[frame].score);
		if (score > target.filter_score) {
			target.filter_score = score;
		}
	}
	//target.filter_time = time;
}

void Query_mapper::greedy_stage(size_t idx, Statistics &stat, int cutoff)
{
	Target& target = targets[idx];
	const sequence subject = ref_seqs::get()[target.subject_id];
	const string subject_id(ref_ids::get()[target.subject_id].c_str());
	if (config.log_subject)
		cout << "Subject = " << subject_id << endl;
#ifdef ENABLE_TIMING
	High_res_timer timer;
#endif
	target.filter_score = 0;
	target.ts.sort(Hsp_traits::cmp_diag);
	typedef Map<list<Hsp_traits>::const_iterator, Hsp_traits::Frame> Hsp_map;
	Hsp_map ts(target.ts.begin(), target.ts.end());
	list<Hsp_traits> t_out;
	target.hsps.clear();
	for (Hsp_map::Iterator it = ts.begin(); it.valid(); ++it) {
		const unsigned frame = it.begin()->frame;
		target.filter_score = std::max(target.filter_score, greedy_align(query_seq(frame), profile[frame], query_cb[frame], subject, log_ga, target.hsps, it.begin(), it.end(), t_out, cutoff, frame));
	}
	target.ts.clear();
	target.ts.splice(target.ts.begin(), t_out);
	//stat.inc(Statistics::TIME_GREEDY_EXT, timer.nanoseconds());
#ifdef ENABLE_TIMING
	target.filter_time = (float)timer.get();
#endif
}

void Query_mapper::align_target(size_t idx, Statistics &stat)
{
	typedef float score_t;
	Target& target = targets[idx];
	const size_t n = target.end - target.begin;
	const vector<Seed_hit>::const_iterator hits = seed_hits.begin() + target.begin;
	const sequence subject = ref_seqs::get()[hits[0].subject_];
	if (config.log_subject)
		cout << "Subject = " << ref_ids::get()[target.subject_id].c_str() << endl;

	stat.inc(Statistics::CELLS, query_seq(0).length() * subject.length());

	if (target.filter_score == 0)
		return;
	if (config.ext == Config::more_greedy)
		target.hsps.push_back(Hsp_data(target.filter_score));
	else {
		const int qlen = (int)query_seq(0).length(),
			band_plus = qlen <= 50 ? 0 : 16;
		target.hsps.clear();
		for (list<Hsp_traits>::const_iterator i = target.ts.begin(); i != target.ts.end(); ++i) {
			if (log_ga) {
				cout << "i_begin=" << i->query_range.begin_ << " j_begin=" << i->subject_range.begin_ << " d_min=" << i->d_min << " d_max=" << i->d_max << endl;
			}
			target.hsps.push_back(Hsp_data());
			target.hsps.back().frame = i->frame;
			banded_sw(query_seq(i->frame), subject, i->d_min - band_plus, i->d_max + band_plus + 1, 0, (int)subject.length(), target.hsps.back());
			if (config.comp_based_stats) {
				const int score = (int)target.hsps.back().score + query_cb[i->frame](target.hsps.back());
				target.hsps.back().score = (unsigned)std::max(0, score);
			}
		}
	}

	if (!target.hsps.empty())
		stat.inc(Statistics::OUT_HITS);

	for (list<Hsp_data>::iterator i = target.hsps.begin(); i != target.hsps.end(); ++i)
		for (list<Hsp_data>::iterator j = target.hsps.begin(); j != target.hsps.end();)
			if (j != i && j->is_weakly_enveloped(*i)) {
				stat.inc(Statistics::ERASED_HITS);
				j = target.hsps.erase(j);
			}
			else
				++j;

	//const float time = (float)timer.getElapsedTimeInMicroSec() + target.filter_time;

	for (list<Hsp_data>::iterator i = target.hsps.begin(); i != target.hsps.end(); ++i) {
		i->time = target.filter_time;
		i->set_source_range(i->frame, source_query_len);
	}

	target.hsps.sort();
	if(target.hsps.size() > 0)
		target.filter_score = target.hsps.front().score;

	target.ts.clear();
	for (list<Hsp_data>::iterator i = target.hsps.begin(); i != target.hsps.end(); ++i)
		target.ts.push_back(Hsp_traits(i->query_source_range));
	
	if (config.use_smith_waterman && !target.hsps.empty()) {
		int score;
		for (unsigned f = 0; f < align_mode.query_contexts; ++f) {
			needleman_wunsch(query_seq(f), subject, score, Local(), int());
			target.hsps.front().sw_score = std::max((unsigned)score, target.hsps.front().sw_score);
		}
		stat.inc(Statistics::SQUARED_ERROR, (stat_type)pow((int)target.hsps.front().sw_score - (int)target.hsps.front().score, 2));
	}
}

void Query_mapper::align_targets(Statistics &stat)
{
	const size_t n = targets.size();
	vector<sequence> seqs(n);
	for (size_t i = 0; i < n; ++i) {
		seqs[i] = ref_seqs::get()[targets[i].subject_id];
	}
	vector<int> scores(n);
	swipe(query_seq(0), seqs.begin(), seqs.end(), scores.begin());
	for (size_t i = 0; i < n; ++i)
		targets[i].hsps.push_back(Hsp_data(scores[i]));
}