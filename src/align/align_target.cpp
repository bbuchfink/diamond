/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include "query_mapper.h"
#include "../data/reference.h"
#include "../dp/floating_sw.h"
#include "../util/map.h"

using std::list;

void Query_mapper::get_prefilter_score(size_t idx)
{
	static const int max_dist = 64;
	static const bool logging = false;

	Target& target = targets[idx];

	if (config.greedy) {
		//cout << "subject=" << ref_ids::get()[(seed_hits.begin() + target.begin)->subject_].c_str() << endl;
		std::sort(seed_hits.begin() + target.begin, seed_hits.begin() + target.end, Seed_hit::compare_diag);
		typedef Map<vector<Seed_hit>::const_iterator, Seed_hit::Frame> Hit_map;		
		Hit_map hit_map(seed_hits.begin() + target.begin, seed_hits.begin() + target.end);
		Hsp_data hsp;
		Hsp_traits traits;
		const sequence subject = ref_seqs::get()[(seed_hits.begin() + target.begin)->subject_];
		for (Hit_map::Iterator it = hit_map.begin(); it.valid(); ++it) {
			const unsigned frame = it.begin()->frame_;
			greedy_align(query_seq(frame), profile[frame], query_cb[frame], subject, it.begin(), it.end(), logging, hsp, traits);
			if (hsp.score > target.filter_score) {
				target.filter_score = hsp.score;
				target.traits = traits;
				target.filter_frame = frame;
				target.filter_i = hsp.query_range.begin_;
				target.filter_j = hsp.subject_range.begin_;
			}
		}
		return;
	}

	const size_t n = target.end - target.begin;
	vector<Seed_hit>::iterator hits = seed_hits.begin() + target.begin;
	std::sort(seed_hits.begin() + target.begin, seed_hits.begin() + target.end, Seed_hit::compare_pos);
	
	int max_score = 0;
	for (unsigned node = 0; node < n; ++node) {
		Seed_hit& d = hits[node];
		if (d.ungapped.len == 0)
			continue;
		for (int k = node - 1; k >= 0; --k) {
			const Seed_hit &e = hits[k];
			if (e.ungapped.len == 0)
				continue;
			if (d.ungapped.j - e.ungapped.subject_last() < max_dist) {
				if (abs(d.ungapped.i - e.ungapped.query_last()) >= max_dist)
					continue;
				const int shift = d.ungapped.diag() - e.ungapped.diag();
				int gap_score = -score_matrix.gap_open() - abs(shift)*score_matrix.gap_extend();
				const int space = shift > 0 ? d.ungapped.j - e.ungapped.subject_last() : d.ungapped.i - e.ungapped.query_last();
				int prefix_score;
				if (space <= 0)
					prefix_score = std::max(e.prefix_score - (e.ungapped.score - e.ungapped.partial_score(abs(space))) + d.ungapped.score, e.prefix_score + d.ungapped.partial_score(abs(space))) + gap_score;
				else
					prefix_score = e.prefix_score + d.ungapped.score + gap_score;

				d.prefix_score = std::max(d.prefix_score, (unsigned)prefix_score);
			}
			else
				break;
		}
		max_score = std::max(max_score, (int)d.prefix_score);
	}
	target.filter_score = max_score;
}

bool is_contained(const vector<Seed_hit>::const_iterator &hits, size_t i)
{
	for (size_t j = 0; j < i; ++j)
		if (hits[i].frame_ == hits[j].frame_ && hits[i].ungapped.is_enveloped(hits[j].ungapped))
			return true;
	return false;
}

bool is_contained(const list<Hsp_data> &hsps, const Seed_hit &hit)
{
	for (list<Hsp_data>::const_iterator i = hsps.begin(); i != hsps.end(); ++i)
		if (hit.frame_ == i->frame && i->pass_through(hit.ungapped))
			return true;
	return false;
}

pair<int, int> get_diag_range(vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end, unsigned frame)
{
	int d_min = std::numeric_limits<int>::max(), d_max = std::numeric_limits<int>::min();
	for (vector<Seed_hit>::const_iterator i = begin; i < end; ++i) {
		if (i->frame_ == frame) {
			const int d = i->diagonal();
			d_min = std::min(d_min, d);
			d_max = std::max(d_max, d);
		}
	}
	return pair<int, int>(d_min, d_max);
}

void Query_mapper::align_target(size_t idx, Statistics &stat)
{
	static const bool logging = false;
	static const int band_plus = 3;
	typedef float score_t;
	Target& target = targets[idx];
	const size_t n = target.end - target.begin,
		max_len = query_seq(0).length() + 100 * query_seqs::get().avg_len();
	size_t aligned_len = 0;
	const vector<Seed_hit>::const_iterator hits = seed_hits.begin() + target.begin;
	const sequence subject = ref_seqs::get()[hits[0].subject_];
	//cout << '>' << query_ids::get()[query_id].c_str() << endl << '>' << ref_ids::get()[hits[0].subject_].c_str() << endl;

	unsigned frame_mask = (1 << align_mode.query_contexts) - 1;
	Timer timer;
	timer.start();

	if (!config.greedy) {

		std::sort(seed_hits.begin() + target.begin, seed_hits.begin() + target.end);

		for (size_t i = 0; i < n; ++i) {
			const unsigned frame = hits[i].frame_;
			if ((frame_mask & (1u << frame)) == 0)
				continue;
			if (!is_contained(hits, i) && !is_contained(target.hsps, hits[i])) {
				target.hsps.push_back(Hsp_data());
				target.hsps.back().frame = frame;
				uint64_t cell_updates;

				if (false && config.comp_based_stats == 1)
					floating_sw(&query_seq(frame)[hits[i].query_pos_],
						&subject[hits[i].subject_pos_],
						target.hsps.back(),
						config.read_padding(query_seq(frame).length()),
						(score_t)score_matrix.rawscore(config.gapped_xdrop),
						(score_t)(score_matrix.gap_open() + score_matrix.gap_extend()),
						(score_t)score_matrix.gap_extend(),
						cell_updates,
						hits[i].query_pos_,
						hits[i].subject_pos_,
						query_cb[frame],
						Traceback(),
						score_t());
				else
					floating_sw(&query_seq(frame)[hits[i].query_pos_],
						&subject[hits[i].subject_pos_],
						target.hsps.back(),
						config.read_padding(query_seq(frame).length()),
						score_matrix.rawscore(config.gapped_xdrop),
						score_matrix.gap_open() + score_matrix.gap_extend(),
						score_matrix.gap_extend(),
						cell_updates,
						hits[i].query_pos_,
						hits[i].subject_pos_,
						No_score_correction(),
						Traceback(),
						int());

				if (config.comp_based_stats) {
					const int score = (int)target.hsps.back().score + query_cb[frame](target.hsps.back());
					target.hsps.back().score = (unsigned)std::max(0, score);
				}

				stat.inc(Statistics::OUT_HITS);
				if (i > 0)
					stat.inc(Statistics::SECONDARY_HITS);
				aligned_len += target.hsps.back().length;
				if (aligned_len > max_len)
					break;
			}
			else
				stat.inc(Statistics::DUPLICATES);
		}
	}
	else {
		if (target.filter_score == 0)
			return;
		const unsigned frame = target.filter_frame;
		const int d = target.filter_i - target.filter_j,
			band = std::max(d - target.traits.d_min, target.traits.d_max - d) + band_plus;
		target.hsps.push_back(Hsp_data());
		target.hsps.back().frame = frame;
		uint64_t cell_updates;
		floating_sw(&query_seq(frame)[target.filter_i],
			&subject[target.filter_j],
			target.hsps.back(),
			band,
			score_matrix.rawscore(config.gapped_xdrop),
			score_matrix.gap_open() + score_matrix.gap_extend(),
			score_matrix.gap_extend(),
			cell_updates,
			target.filter_i,
			target.filter_j,
			No_score_correction(),
			Traceback(),
			int());

		if (config.comp_based_stats) {
			const int score = (int)target.hsps.back().score + query_cb[frame](target.hsps.back());
			target.hsps.back().score = (unsigned)std::max(0, score);
		}

		stat.inc(Statistics::OUT_HITS);
	}

	for (list<Hsp_data>::iterator i = target.hsps.begin(); i != target.hsps.end(); ++i)
		for (list<Hsp_data>::iterator j = target.hsps.begin(); j != target.hsps.end();)
			if (j != i && j->is_weakly_enveloped(*i)) {
				stat.inc(Statistics::ERASED_HITS);
				j = target.hsps.erase(j);
			}
			else
				++j;

	const float time = (float)timer.getElapsedTimeInMicroSec();

	for (list<Hsp_data>::iterator i = target.hsps.begin(); i != target.hsps.end(); ++i) {
		i->time = time;
		i->set_source_range(i->frame, source_query_len);
	}

	target.hsps.sort();
	if(target.hsps.size() > 0)
		target.filter_score = target.hsps.front().score;
}