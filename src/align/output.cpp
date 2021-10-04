/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
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
#include "../output/daa/daa_write.h"
#include "../util/sequence/sequence.h"

using std::vector;

namespace Extension {

TextBuffer* generate_output(vector<Match> &targets, size_t query_block_id, Statistics &stat, const Search::Config& cfg)
{
	const SequenceSet& query_seqs = cfg.query->seqs(), &ref_seqs = cfg.target->seqs();
	TextBuffer* out = new TextBuffer;
	std::unique_ptr<Output_format> f(output_format->clone());
	size_t seek_pos = 0;
	unsigned n_hsp = 0, hit_hsps = 0;
	TranslatedSequence query = query_seqs.translated_seq(align_mode.query_translated ? cfg.query->source_seqs()[query_block_id] : query_seqs[query_block_id], query_block_id * align_mode.query_contexts);
	const unsigned query_len = (unsigned)query.index(0).length();
	const char *query_title = cfg.query->ids()[query_block_id];
	const double query_self_aln_score = cfg.query->has_self_aln() ? cfg.query->self_aln_score(query_block_id) : 0.0;
	const bool aligned = !targets.empty();

	if (cfg.iterated()) {
		if (aligned) seek_pos = IntermediateRecord::write_query_intro(*out, query_block_id);
	}
	else if (*f == Output_format::daa) {
		if (aligned) seek_pos = write_daa_query_record(*out, query_title, query.source());
	}
	else if (aligned || config.report_unaligned)
		f->print_query_intro(query_block_id, query_title, query.source().length(), *out, !aligned, cfg);
			
	for (size_t i = 0; i < targets.size(); ++i) {

		if (targets[i].hsp.empty())
			throw std::runtime_error("generate_output: target with no hsps.");

		const size_t subject_id = targets[i].target_block_id;
		const unsigned database_id = cfg.target->block_id2oid(subject_id);
		const unsigned subject_len = (unsigned)ref_seqs[subject_id].length();
		const double target_self_aln_score = cfg.target->has_self_aln() ? cfg.target->self_aln_score(subject_id) : 0.0;

		hit_hsps = 0;
		for (Hsp &hsp : targets[i].hsp) {
			if (*f == Output_format::daa)
				write_daa_record(*out, hsp, cfg.target->dict_id(current_ref_block, subject_id, *cfg.db));
			else if(cfg.iterated())
				IntermediateRecord::write(*out, hsp, query_block_id, cfg.target->dict_id(current_ref_block, subject_id, *cfg.db), database_id);
			else
				f->print_match(HspContext(hsp,
					query_block_id,
					query,
					query_title,
					database_id,
					subject_len,
					(cfg.target->has_ids() ? cfg.target->ids()[subject_id] : cfg.db->seqid(database_id)).c_str(),
					i,
					hit_hsps,
					cfg.target->unmasked_seqs().empty() ? Sequence() : cfg.target->unmasked_seqs()[subject_id],
					targets[i].ungapped_score,
					query_self_aln_score,
					target_self_aln_score), cfg, *out);
			
			++n_hsp;
			++hit_hsps;
		}
	}

	if (cfg.iterated()) {
		if (aligned) IntermediateRecord::finish_query(*out, seek_pos);
	}
	else {
		stat.inc(Statistics::MATCHES, n_hsp);
		stat.inc(Statistics::PAIRWISE, targets.size());
		if (aligned)
			stat.inc(Statistics::ALIGNED);
		if (*f == Output_format::daa) {
			if (aligned) finish_daa_query_record(*out, seek_pos);
		}
		else if (aligned || config.report_unaligned)
			f->print_query_epilog(*out, query_title, targets.empty(), cfg);
	}

	return out;
}

TextBuffer* generate_intermediate_output(const vector<Match> &targets, size_t query_block_id, const Search::Config& cfg)
{
	TextBuffer* out = new TextBuffer;
	if (targets.empty())
		return out;
	size_t seek_pos = 0;
	seek_pos = IntermediateRecord::write_query_intro(*out, query_block_id);
	const Block& target = *cfg.target;
	
	for (size_t i = 0; i < targets.size(); ++i) {

		const size_t block_id = targets[i].target_block_id;
		const size_t dict_id = target.dict_id(current_ref_block, block_id, *cfg.db);

		for (const Hsp &hsp : targets[i].hsp)
			IntermediateRecord::write(*out, hsp, query_block_id, dict_id, target.block_id2oid(block_id));
		/*if (config.global_ranking_targets > 0)
			IntermediateRecord::write(*out, subject_id, targets[i].ungapped_score, cfg);*/
	}

	IntermediateRecord::finish_query(*out, seek_pos);
	return out;
}

}
