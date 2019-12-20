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

#include "extend.h"
#include "target.h"
#include "../data/queries.h"

using std::vector;

namespace Extension {

bool generate_output(const vector<Match> &targets, size_t query_id, TextBuffer &buffer, Statistics &stat)
{
	unsigned n_hsp = 0, n_target_seq = 0, hit_hsps = 0;
	TranslatedSequence query = query_seqs::get().translated_seq(query_source_seqs::get()[query_id], query_id*align_mode.query_contexts);

	/*const unsigned query_len = (unsigned)query_seq(0).length();
	size_t seek_pos = 0;
	const char *query_title = query_ids::get()[query_id].c_str();
	unique_ptr<Output_format> f(output_format->clone());

	for (size_t i = 0; i < targets.size(); ++i) {

		if ((config.min_bit_score == 0 && score_matrix.evalue(targets[i].filter_score, query_len) > config.max_evalue)
			|| score_matrix.bitscore(targets[i].filter_score) < config.min_bit_score)
			break;

		const size_t subject_id = targets[i].subject_block_id;
		const unsigned database_id = ReferenceDictionary::get().block_to_database_id(subject_id);
		const unsigned subject_len = (unsigned)ref_seqs::get()[subject_id].length();
		const char *ref_title = ref_ids::get()[subject_id].c_str();
		targets[i].apply_filters(source_query_len, subject_len, query_title, ref_title);
		if (targets[i].hsps.size() == 0)
			continue;

		const int c = target_culling->cull(targets[i]);
		if (c == TargetCulling::NEXT)
			continue;
		else if (c == TargetCulling::FINISHED)
			break;

		if (targets[i].outranked)
			stat.inc(Statistics::OUTRANKED_HITS);
		target_culling->add(targets[i]);

		hit_hsps = 0;
		for (list<Hsp>::iterator j = targets[i].hsps.begin(); j != targets[i].hsps.end(); ++j) {
			if (config.max_hsps > 0 && hit_hsps >= config.max_hsps)
				break;

			if (blocked_processing) {
				if (n_hsp == 0)
					seek_pos = IntermediateRecord::write_query_intro(buffer, query_id);
				IntermediateRecord::write(buffer, *j, query_id, subject_id);
			}
			else {
				if (n_hsp == 0) {
					if (*f == Output_format::daa)
						seek_pos = write_daa_query_record(buffer, query_title, align_mode.query_translated ? query_source_seqs::get()[query_id] : query_seqs::get()[query_id]);
					else
						f->print_query_intro(query_id, query_title, source_query_len, buffer, false);
				}
				if (*f == Output_format::daa)
					write_daa_record(buffer, *j, subject_id);
				else
					f->print_match(Hsp_context(*j,
						query_id,
						translated_query,
						query_title,
						subject_id,
						database_id,
						ref_title,
						subject_len,
						n_target_seq,
						hit_hsps,
						ref_seqs::get()[subject_id]), metadata, buffer);
			}

			++n_hsp;
			++hit_hsps;
		}
		++n_target_seq;
	}

	if (n_hsp > 0) {
		if (!blocked_processing) {
			if (*f == Output_format::daa)
				finish_daa_query_record(buffer, seek_pos);
			else
				f->print_query_epilog(buffer, query_title, false, parameters);
		}
		else
			IntermediateRecord::finish_query(buffer, seek_pos);
	}
	else if (!blocked_processing && *f != Output_format::daa && config.report_unaligned != 0) {
		f->print_query_intro(query_id, query_title, source_query_len, buffer, true);
		f->print_query_epilog(buffer, query_title, true, parameters);
	}

	if (!blocked_processing) {
		stat.inc(Statistics::MATCHES, n_hsp);
		stat.inc(Statistics::PAIRWISE, n_target_seq);
		if (n_hsp > 0)
			stat.inc(Statistics::ALIGNED);
	}*/

	return n_hsp > 0;
}

}