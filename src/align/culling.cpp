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

#include <algorithm>
#include "target.h"
#include "../basic/config.h"
#include "../data/reference.h"
#include "culling.h"

using std::vector;
using std::list;

namespace Extension {

static void max_hsp_culling(list<Hsp>& hsps) {
	if (config.max_hsps > 0 && hsps.size() > config.max_hsps) {
		list<Hsp>::iterator i = hsps.begin();
		for (unsigned n = 0; n < config.max_hsps; ++n)
			++i;
		hsps.erase(i, hsps.end());
	}
}

static void inner_culling(list<Hsp>& hsps, int source_query_len) {
	for (Hsp& h : hsps)
		h.query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(h.query_range.begin_, Frame(h.frame)), TranslatedPosition(h.query_range.end_, Frame(h.frame)), source_query_len);
	hsps.sort();
	const double overlap = config.inner_culling_overlap / 100.0;
	for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end();) {
		if (i->is_enveloped_by(hsps.begin(), i, overlap))
			i = hsps.erase(i);
		else
			++i;
	}
	if (config.max_hsps > 0)
		max_hsp_culling(hsps);
}

void Target::inner_culling(int source_query_len) {
	list<Hsp> hsps;
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame)
		hsps.splice(hsps.end(), hsp[frame]);
	Extension::inner_culling(hsps, source_query_len);
	while (!hsps.empty()) {
		auto& l = hsp[hsps.front().frame];
		l.splice(l.end(), hsps, hsps.begin());
	}
}

void Match::inner_culling(int source_query_len)
{
	Extension::inner_culling(hsp, source_query_len);
	if (!hsp.empty())
		filter_score = hsp.front().score;
}

void Match::max_hsp_culling() {
	Extension::max_hsp_culling(hsp);
}

bool append_hits(vector<Target>& targets, vector<Target>::const_iterator begin, vector<Target>::const_iterator end, int low_score, size_t previous_count, int source_query_len, const char* query_title, const sequence& query_seq) {
	bool append = false;

	if (config.toppercent != 100.0 || targets.size() >= config.max_alignments) {
		culling(targets, source_query_len, query_title, query_seq);
	}
	else
		append = true;

	int max_score = 0;
	for (auto i = begin; i < end; ++i)
		max_score = std::max(max_score, i->filter_score);

	if (targets.empty())
		append = true;
	else if (config.toppercent == 100.0 && max_score >= targets.back().filter_score)
		append = true;
	else if (config.toppercent != 100.0 && max_score >= top_cutoff_score(targets.front().filter_score))
		append = true;

	/*if (max_score <= low_score && previous_count >= config.max_alignments && config.toppercent == 100.0)
		return false;*/

	/*if (config.toppercent == 100.0 && (config.min_id > 0 || config.query_cover > 0 || config.subject_cover > 0 || config.no_self_hits))
		return true;*/	

	if(append)
		targets.insert(targets.end(), begin, end);

	return append;
}

/*void score_only_culling(vector<Target> &targets) {
	const size_t n = targets.size();	

	std::sort(targets.begin(), targets.end());
	if (targets.empty() || targets.front().filter_score == 0) {
		targets.clear();
		return;
	}

	if (config.toppercent == 100.0 && (config.min_id > 0 || config.query_cover > 0 || config.subject_cover > 0 || config.no_self_hits))
		return;

	vector<Target>::iterator i = targets.begin();
	if (config.toppercent < 100.0) {
		const int cutoff = std::max(int((1.0 - config.toppercent / 100.0)*targets.front().filter_score), 1);
		while (i < targets.end() && i->filter_score >= cutoff)
			++i;
	}
	else {
		i += std::min((size_t)config.max_alignments, targets.size());
		while (--i > targets.begin() && i->filter_score == 0);
		++i;
	}
	targets.erase(i, targets.end());
}*/

bool filter_hsp(const Hsp& hsp, int source_query_len, const char *query_title, int subject_len, const char* subject_title, const sequence& query_seq, const sequence& subject_seq) {
	return 	hsp.id_percent() < config.min_id
		|| hsp.query_cover_percent(source_query_len) < config.query_cover
		|| hsp.subject_cover_percent(subject_len) < config.subject_cover
		|| (config.no_self_hits
			&& query_seq == subject_seq
			&& strcmp(query_title, subject_title) == 0);
}

void Target::apply_filters(int source_query_len, const char *query_title, const sequence& query_seq)
{
	const char *title = ref_ids::get()[block_id];
	const sequence seq = ref_seqs::get()[block_id];
	const int len = seq.length();
	filter_score = 0;
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		for (list<Hsp>::iterator i = hsp[frame].begin(); i != hsp[frame].end();) {
			if (filter_hsp(*i, source_query_len, query_title, len, title, query_seq, seq))
				i = hsp[frame].erase(i);
			else {
				++i;
				filter_score = std::max(filter_score, i->score);
			}
		}
	}
}

void Match::apply_filters(int source_query_len, const char *query_title, const sequence& query_seq)
{
	const char *title = ref_ids::get()[target_block_id];
	const sequence seq = ref_seqs::get()[target_block_id];
	const int len = seq.length();
	for (list<Hsp>::iterator i = hsp.begin(); i != hsp.end();) {
		if (filter_hsp(*i, source_query_len, query_title, len, title, query_seq, seq))
			i = hsp.erase(i);
		else
			++i;
	}
	filter_score = hsp.empty() ? 0 : hsp.front().score;
}

/*void culling(vector<Match> &targets, int source_query_len, const char *query_title) {
	for (Match &match : targets)
		match.apply_filters(source_query_len, query_title);
	std::sort(targets.begin(), targets.end());
	if (targets.empty() || targets.front().filter_score == 0) {
		targets.clear();
		return;
	}
	vector<Match>::iterator i = targets.begin();
	if (config.toppercent < 100.0) {
		const int cutoff = std::max(int((1.0 - config.toppercent / 100.0)*targets.front().filter_score), 1);
		while (i < targets.end() && i->filter_score >= cutoff)
			++i;
	}
	else {
		i += std::min((size_t)config.max_alignments, targets.size());
		while (--i > targets.begin() && i->filter_score == 0);
		++i;
	}
	targets.erase(i, targets.end());
}*/

}