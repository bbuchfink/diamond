/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

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

#include <algorithm>
#include <thread>
#include <atomic>
#include "../align.h"
#include "../../dp/dp.h"
#include "../../util/geo/interval_partition.h"
#include "../../util/simd.h"

using std::thread;
using std::list;
using std::max;
using std::atomic;
using std::endl;
using std::vector;

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

	Interval ungapped_query_range(int query_dna_len) const
	{
		const Frame f = Frame(top_hit.frame_);
		const int i0 = std::max((int)top_hit.query_pos_ - (int)top_hit.subject_pos_, 0),
			i1 = std::min((int)top_hit.query_pos_ + (int)subject.length() - (int)top_hit.subject_pos_, f.length(query_dna_len));
		return TranslatedPosition::absolute_interval(TranslatedPosition(i0, f), TranslatedPosition(i1, f), query_dna_len);
	}

	void add_strand(QueryMapper &mapper, vector<DpTarget> &v, vector<Seed_hit>::iterator begin, vector<Seed_hit>::iterator end, int target_idx)
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
				v.emplace_back(subject, subject.length(), d0, d1, Interval(), 0, target_idx, 0); // set cols here?
				d0 = std::max(i->diagonal() - band, d_min);
				d1 = std::min(i->diagonal() + band, d_max);
			}
		}
		v.emplace_back(subject, subject.length(), d0, d1, Interval(), 0, target_idx, 0);
	}

	void add(QueryMapper &mapper, vector<DpTarget> &vf, vector<DpTarget> &vr, int target_idx)
	{
		vector<Seed_hit>::iterator hits = mapper.seed_hits.begin() + begin, hits_end = mapper.seed_hits.begin() + end;
		Strand strand = top_hit.strand();
		if(mapper.target_parallel)
			std::stable_sort(hits, hits_end, Seed_hit::compare_diag_strand);
		else
			std::stable_sort(hits, hits_end, Seed_hit::compare_diag_strand2);

		const auto it = find_if(hits, hits_end, [](const Seed_hit &x) { return x.strand() == REVERSE; });
		if(strand == FORWARD || mapper.target_parallel)
			add_strand(mapper, vf, hits, it, target_idx);
		if (strand == REVERSE || mapper.target_parallel)
			add_strand(mapper, vr, it, hits_end, target_idx);

		//const int d = hits[0].diagonal();
		//const DpTarget t(subject, std::max(d - 32, -int(subject.length() - 1)), std::min(d + 32, int(mapper.query_seq(0).length() - 1)), &hsps, subject_id);		
	}

	void set_filter_score()
	{
		filter_score = 0;
		filter_evalue = DBL_MAX;
		for (list<Hsp>::const_iterator i = hsps.begin(); i != hsps.end(); ++i) {
			filter_score = std::max(filter_score, (int)i->score);
			filter_evalue = std::min(filter_evalue, i->evalue);
		}
	}

	void reset()
	{
		hsps.clear();
	}

	void finish(QueryMapper &mapper)
	{
		inner_culling();
		if (config.frame_shift)
			return;
		for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end(); ++i)
			i->query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(i->query_range.begin_, Frame(i->frame)), TranslatedPosition(i->query_range.end_, Frame(i->frame)), mapper.source_query_len);
	}

	bool is_outranked(const IntervalPartition &ip, int source_query_len, double rr) const
	{
		const Interval r = ungapped_query_range(source_query_len);
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

void Pipeline::range_ranking(const int64_t max_target_seqs)
{
	const double rr = config.rank_ratio == -1 ? 0.4 : config.rank_ratio;
	std::stable_sort(targets.begin(), targets.end(), Target::compare_score);
	IntervalPartition ip(max_target_seqs);
	for (PtrVector< ::Target>::iterator i = targets.begin(); i < targets.end();) {
		Target* t = ((Target*)*i);
		if (t->is_outranked(ip, source_query_len, rr)) {
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
	for (int64_t i = 0; i < n_targets(); ++i)
		target(i).add(*this, vf, vr, (int)i);
	list<Hsp> hsp;
	hsp = banded_3frame_swipe(translated_query, FORWARD, vf.begin(), vf.end(), this->dp_stat, score_only, target_parallel);
	hsp.splice(hsp.end(), banded_3frame_swipe(translated_query, REVERSE, vr.begin(), vr.end(), this->dp_stat, score_only, target_parallel));
	
	while (!hsp.empty()) {
		list<Hsp> &l = target(hsp.begin()->swipe_target).hsps;
		l.splice(l.end(), hsp, hsp.begin());
	}
}

void build_ranking_worker(PtrVector<::Target>::iterator begin, PtrVector<::Target>::iterator end, atomic<size_t> *next, vector<int32_t> *intervals) {
	PtrVector<::Target>::iterator it;
	while ((it = begin + next->fetch_add(64)) < end) {
		auto e = min(it + 64, end);
		for (; it < e; ++it) {
			(*it)->add_ranges(*intervals);
		}
	}
}

void Pipeline::run(Statistics &stat, const Search::Config& cfg)
{
	task_timer timer("Init banded swipe pipeline", target_parallel ? 3 : UINT_MAX);
	Config::set_option(config.padding, 32);
	if (n_targets() == 0)
		return;
	stat.inc(Statistics::TARGET_HITS0, n_targets());

	if (!target_parallel) {
		timer.go("Ungapped stage");
		for (int64_t i = 0; i < n_targets(); ++i)
			target(i).ungapped_stage(*this);
		timer.go("Ranking");
		if (!config.query_range_culling)
			rank_targets(config.rank_ratio == -1 ? 0.4 : config.rank_ratio, config.rank_factor == -1.0 ? 1e3 : config.rank_factor, cfg.max_target_seqs);
		else
			range_ranking(cfg.max_target_seqs);
	}
	else {
		timer.finish();
		log_stream << "Query: " << query_id << "; Seed hits: " << seed_hits.size() << "; Targets: " << n_targets() << endl;
	}

	if (n_targets() > cfg.max_target_seqs || config.toppercent < 100.0) {
		stat.inc(Statistics::TARGET_HITS1, n_targets());
		timer.go("Swipe (score only)");
		run_swipe(true);

		if (target_parallel) {
			timer.go("Building score ranking intervals");
			vector<vector<int32_t>> intervals(config.threads_);
			const size_t interval_count = (source_query_len + ::Target::INTERVAL - 1) / ::Target::INTERVAL;
			for (vector<int32_t> &v : intervals)
				v.resize(interval_count);
			vector<thread> threads;
			atomic<size_t> next(0);
			for (int i = 0; i < config.threads_; ++i)
				threads.emplace_back(build_ranking_worker, targets.begin(), targets.end(), &next, &intervals[i]);
			for (auto &t : threads)
				t.join();

			timer.go("Merging score ranking intervals");
			for (auto it = intervals.begin() + 1; it < intervals.end(); ++it) {
				vector<int32_t>::iterator i = intervals[0].begin(), i1 = intervals[0].end(), j = it->begin();
				for (; i < i1; ++i, ++j)
					*i = max(*i, *j);
			}

			timer.go("Finding outranked targets");
			for (auto i = targets.begin(); i < targets.end(); ++i)
				if ((*i)->is_outranked(intervals[0], 1.0 - config.toppercent / 100)) {
					delete *i;
					*i = nullptr;
				}

			timer.go("Removing outranked targets");
			vector<::Target*> &v(*static_cast<vector<::Target*>*>(&targets));
			auto it = remove_if(v.begin(), v.end(), [](::Target *t) { return t == nullptr; });
			v.erase(it, v.end());
			timer.finish();
			log_stream << "Targets after score-only ranking: " << targets.size() << endl;
		}
		else {
			timer.go("Score only culling");
			for (int64_t i = 0; i < n_targets(); ++i)
				target(i).set_filter_score();
			score_only_culling(cfg.max_target_seqs);
		}
	}

	timer.go("Swipe (traceback)");
	stat.inc(Statistics::TARGET_HITS2, n_targets());
	for (int64_t i = 0; i < n_targets(); ++i)
		target(i).reset();
	run_swipe(false);

	timer.go("Inner culling");
	for (int64_t i = 0; i < n_targets(); ++i)
		target(i).finish(*this);
}

}}
