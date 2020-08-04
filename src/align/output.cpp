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

#include <memory>
#include "extend.h"
#include "target.h"
#include "../data/queries.h"
#include "../output/output_format.h"
#include "../data/ref_dictionary.h"
#include "../output/daa_write.h"

using std::vector;

namespace Extension {

TextBuffer* generate_output(vector<Match> &targets, size_t query_block_id, Statistics &stat, const Metadata &metadata, const Parameters &parameters)
{
	TextBuffer* out = new TextBuffer;
	std::unique_ptr<Output_format> f(output_format->clone());
	size_t seek_pos = 0;
	unsigned n_hsp = 0, hit_hsps = 0;
	TranslatedSequence query = query_seqs::get().translated_seq(align_mode.query_translated ? query_source_seqs::get()[query_block_id] : query_seqs::get()[query_block_id], query_block_id*align_mode.query_contexts);
	const unsigned query_len = (unsigned)query.index(0).length();
	const char *query_title = query_ids::get()[query_block_id];
	const bool aligned = !targets.empty();

	if (blocked_processing) {
		if(aligned) seek_pos = IntermediateRecord::write_query_intro(*out, query_block_id);
	} else {
		if (*f == Output_format::daa) {
			if (aligned) seek_pos = write_daa_query_record(*out, query_title, query.source());
		} else if(aligned || config.report_unaligned)
			f->print_query_intro(query_block_id, query_title, query.source().length(), *out, !aligned);
	}
	
	for (size_t i = 0; i < targets.size(); ++i) {

		const size_t subject_id = targets[i].target_block_id;
		const unsigned database_id = ReferenceDictionary::get().block_to_database_id(subject_id);
		const unsigned subject_len = (unsigned)ref_seqs::get()[subject_id].length();
		const char *ref_title = ref_ids::get()[subject_id];

		if (targets[i].outranked)
			stat.inc(Statistics::OUTRANKED_HITS);

		hit_hsps = 0;
		for (Hsp &hsp : targets[i].hsp) {
			if (blocked_processing) {
				IntermediateRecord::write(*out, hsp, query_block_id, subject_id);
			}
			else {
				if (*f == Output_format::daa)
					write_daa_record(*out, hsp, subject_id);
				else
					f->print_match(Hsp_context(hsp,
						query_block_id,
						query,
						query_title,
						subject_id,
						database_id,
						ref_title,
						subject_len,
						i,
						hit_hsps,
						ref_seqs::get()[subject_id],
						targets[i].ungapped_score), metadata, *out);
			}

			++n_hsp;
			++hit_hsps;
		}
	}

	if (!blocked_processing) {
		if (*f == Output_format::daa) {
			if(aligned) finish_daa_query_record(*out, seek_pos);
		} else if(aligned || config.report_unaligned)
			f->print_query_epilog(*out, query_title, targets.empty(), parameters);
	}
	else if(aligned)
		IntermediateRecord::finish_query(*out, seek_pos);

	if (!blocked_processing) {
		stat.inc(Statistics::MATCHES, n_hsp);
		stat.inc(Statistics::PAIRWISE, targets.size());
		if (aligned)
			stat.inc(Statistics::ALIGNED);
	}
	return out;
}

}