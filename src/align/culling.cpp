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

using std::vector;
using std::list;

namespace Extension {

static void max_hsp_culling(list<Hsp>& hsps) {
	if (config.max_hsps > 0 && hsps.size() > config.max_hsps)
		hsps.resize(config.max_hsps);
}

static void inner_culling(list<Hsp>& hsps) {
	if (hsps.size() <= 1)
		return;
	hsps.sort();
	if (config.max_hsps == 1) {
		hsps.resize(1);
		return;
	}
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

void Target::inner_culling() {
	if (config.max_hsps == 1) {
		for (int i = 0; i < MAX_CONTEXT; ++i)
			if (i == best_context) {
				hsp[i].sort();
				hsp[i].resize(1);
			}
			else
				hsp[i].clear();			
		return;
	}
	list<Hsp> hsps;
	for (int frame = 0; frame < align_mode.query_contexts; ++frame)
		hsps.splice(hsps.end(), hsp[frame]);
	Extension::inner_culling(hsps);
	while (!hsps.empty()) {
		auto& l = hsp[hsps.front().frame];
		l.splice(l.end(), hsps, hsps.begin());
	}
}

void Match::inner_culling()
{
	Extension::inner_culling(hsp);
	if (!hsp.empty()) {
		filter_evalue = hsp.front().evalue;
		filter_score = hsp.front().score;
	}
}

void Match::max_hsp_culling() {
	Extension::max_hsp_culling(hsp);
}

bool append_hits(vector<Target>& targets, vector<Target>::const_iterator begin, vector<Target>::const_iterator end, size_t chunk_size, int source_query_len, const char* query_title, const Sequence& query_seq, const Block& target_block) {
	bool append = false;

	if (config.toppercent != 100.0 || targets.size() >= config.max_alignments) {
		culling(targets, source_query_len, query_title, query_seq, 0, target_block);
	}
	else
		append = true;

	int max_score = 0;
	double min_evalue = DBL_MAX;
	for (auto i = begin; i < end; ++i) {
		max_score = std::max(max_score, i->filter_score);
		min_evalue = std::min(min_evalue, i->filter_evalue);
	}

	if (targets.empty())
		append = true;
	else if (config.toppercent == 100.0 && min_evalue <= targets.back().filter_evalue)
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

bool filter_hsp(const Hsp& hsp, int source_query_len, const char *query_title, int subject_len, const char* subject_title, const Sequence& query_seq, const Sequence& subject_seq) {
	return 	hsp.id_percent() < config.min_id
		|| hsp.query_cover_percent(source_query_len) < config.query_cover
		|| hsp.subject_cover_percent(subject_len) < config.subject_cover
		|| (config.no_self_hits
			&& query_seq == subject_seq
			&& strcmp(query_title, subject_title) == 0);
}

void Target::apply_filters(int source_query_len, const char *query_title, const Sequence& query_seq, const Block& targets)
{
	const char* title = config.no_self_hits ? targets.ids()[block_id] : nullptr;
	const Sequence seq = targets.seqs()[block_id];
	const int len = seq.length();
	filter_score = 0;
	filter_evalue = DBL_MAX;
	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		for (list<Hsp>::iterator i = hsp[frame].begin(); i != hsp[frame].end();) {
			if (filter_hsp(*i, source_query_len, query_title, len, title, query_seq, seq))
				i = hsp[frame].erase(i);
			else {
				filter_score = std::max(filter_score, i->score);
				filter_evalue = std::min(filter_evalue, i->evalue);
				++i;
			}
		}
	}
}

void Match::apply_filters(int source_query_len, const char *query_title, const Sequence& query_seq, const Block& targets)
{
	const char* title = config.no_self_hits ? targets.ids()[target_block_id] : nullptr;
	const Sequence seq = targets.seqs()[target_block_id];
	const int len = seq.length();
	for (list<Hsp>::iterator i = hsp.begin(); i != hsp.end();) {
		if (filter_hsp(*i, source_query_len, query_title, len, title, query_seq, seq))
			i = hsp.erase(i);
		else
			++i;
	}
	filter_evalue = hsp.empty() ? DBL_MAX : hsp.front().evalue;
	filter_score = hsp.empty() ? 0 : hsp.front().score;
}

void culling(std::vector<Target>& targets, int source_query_len, const char* query_title, const Sequence& query_seq, size_t min_keep, const Block& target_block) {
	if (config.min_id > 0 || config.query_cover > 0 || config.subject_cover > 0 || config.no_self_hits)
		for (Target& match : targets)
			match.apply_filters(source_query_len, query_title, query_seq, target_block);

	std::sort(targets.begin(), targets.end(), config.toppercent < 100.0 ? Target::comp_score : Target::comp_evalue);
	if (targets.empty() || targets.front().filter_evalue == DBL_MAX) {
		targets.clear();
		return;
	}

	std::vector<Target>::iterator i = targets.begin();
	if (config.toppercent < 100.0) {
		size_t n = 0;
		const double cutoff = std::max(top_cutoff_score(score_matrix.bitscore(targets.front().filter_score)), 1.0);
		while (i < targets.end() && (score_matrix.bitscore(i->filter_score) >= cutoff || n < min_keep)) {
			++i;
			++n;
		}
	}
	else {
		i += std::min(std::max(config.max_alignments, min_keep), targets.size());
		while (--i > targets.begin() && i->filter_evalue == DBL_MAX);
		++i;
	}
	targets.erase(i, targets.end());
}

}