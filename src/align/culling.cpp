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

void Match::inner_culling(int source_query_len)
{
	hsp.sort();
	if (!hsp.empty())
		filter_score = hsp.front().score;
	for(Hsp &h : hsp)
		h.query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(h.query_range.begin_, Frame(h.frame)), TranslatedPosition(h.query_range.end_, Frame(h.frame)), source_query_len);
	const double overlap = config.inner_culling_overlap / 100.0;
	for (list<Hsp>::iterator i = hsp.begin(); i != hsp.end();) {
		if (i->is_enveloped_by(hsp.begin(), i, overlap))
			i = hsp.erase(i);
		else
			++i;
	}
	if (config.max_hsps > 0)
		max_hsp_culling();
}

void Match::max_hsp_culling() {
	if (hsp.size() > config.max_hsps) {
		list<Hsp>::iterator i = hsp.begin();
		for (unsigned n = 0; n < config.max_hsps; ++n)
			++i;
		hsp.erase(i, hsp.end());
	}
}

void score_only_culling(vector<Target> &targets) {
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
}

void Match::apply_filters(int source_query_len, const char *query_title)
{
	const char *title = ref_ids::get()[target_block_id];
	const int len = (int)ref_seqs::get()[target_block_id].length();
	for (list<Hsp>::iterator i = hsp.begin(); i != hsp.end();) {
		if (i->id_percent() < config.min_id
			|| i->query_cover_percent(source_query_len) < config.query_cover
			|| i->subject_cover_percent(len) < config.subject_cover
			|| (config.no_self_hits
				&& i->identities == i->length
				&& i->query_source_range.length() == source_query_len
				&& i->subject_range.length() == len
				&& strcmp(query_title, title) == 0))
			i = hsp.erase(i);
		else
			++i;
	}
	filter_score = hsp.empty() ? 0 : hsp.front().score;
}

void culling(vector<Match> &targets, int source_query_len, const char *query_title) {
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
}

}