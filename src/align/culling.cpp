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
	for (list<Hsp>::iterator i = hsp.begin(); i != hsp.end();) {
		if (i->is_enveloped_by(hsp.begin(), i, 0.5))
			i = hsp.erase(i);
		else
			++i;
	}
}

void score_only_culling(vector<Target> &targets) {
	std::stable_sort(targets.begin(), targets.end());

	if (config.toppercent == 100.0 && (config.min_id > 0 || config.query_cover > 0 || config.subject_cover > 0 || config.no_self_hits))
		return;

	vector<Target>::iterator i;
	if (config.toppercent < 100.0) {
		const int cutoff = (1.0 - config.toppercent / 100.0)*targets.front().filter_score;
		while (i < targets.end() && i->filter_score >= cutoff)
			++i;
	}
	else
		i = targets.begin() + std::min((size_t)config.max_alignments, targets.size());
	targets.erase(i, targets.end());
}

void Match::apply_filters(int source_query_len, const char *query_title)
{
	/*for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end();) {
		if (i->id_percent() < config.min_id
			|| i->query_cover_percent(dna_len) < config.query_cover
			|| i->subject_cover_percent(subject_len) < config.subject_cover
			|| (config.no_self_hits
				&& i->identities == i->length
				&& i->query_source_range.length() == (int)dna_len
				&& i->subject_range.length() == (int)subject_len
				&& strcmp(query_title, ref_title) == 0)
			|| (config.filter_locus && !i->subject_range.includes(config.filter_locus)))
			i = hsps.erase(i);
		else
			++i;
	}*/
}

void culling(vector<Match> &targets) {
	std::stable_sort(targets.begin(), targets.end());
}

}