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
#include <limits.h>
#include "extend.h"
#include "../data/queries.h"
#include "../basic/config.h"
#include "../dp/comp_based_stats.h"
#include "target.h"
#include "../data/reference.h"
#include "../data/ref_dictionary.h"

using std::vector;

namespace Extension {

vector<Target> init_targets(sequence *query_seq, Trace_pt_list::iterator begin, Trace_pt_list::iterator end, const Metadata &metadata) {
	vector<Target> targets;
	std::sort(begin, end, hit::cmp_subject);
	size_t target = SIZE_MAX;
	Trace_pt_list::iterator target_begin = begin;
	for (Trace_pt_list::iterator i = begin; i < end; ++i) {
		std::pair<size_t, size_t> l = ref_seqs::data_->local_position(i->subject_);
		if (l.first != target) {
			if (i - target_begin > 0)
				targets.emplace_back(
					target_begin,
					i,
					ref_seqs::data_->position(l.first, 0),
					query_seq,
					ref_seqs::get()[target],
					config.taxon_k ? metadata.taxon_nodes->rank_taxid((*metadata.taxon_list)[ReferenceDictionary::get().block_to_database_id(target)], Rank::species) : set<unsigned>());
			target = l.first;
			target_begin = i;
		}
	}
	return targets;
}

void extend(const Parameters &params, size_t query_id, Trace_pt_list::iterator begin, Trace_pt_list::iterator end, const Metadata &metadata, Statistics &stat) {
	const unsigned contexts = align_mode.query_contexts;
	vector<sequence> query_seq;
	vector<Bias_correction> query_cb;

	if (config.log_query)
		cout << "Query = " << query_ids::get()[query_id].c_str() << endl;

	for (unsigned i = 0; i < contexts; ++i)
		query_seq.push_back(query_seqs::get()[query_id*contexts + i]);

	if (config.comp_based_stats == 1)
		for (unsigned i = 0; i < align_mode.query_contexts; ++i)
			query_cb.emplace_back(query_seq[i]);

	vector<Target> targets = init_targets(query_seq.data(), begin, end, metadata);
	stat.inc(Statistics::TARGET_HITS0, targets.size());

}

}