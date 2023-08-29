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

#include <algorithm>
#include "target.h"
#include "../basic/config.h"
#include "../data/reference.h"
#include "../output/recursive_parser.h"
#include "../output/output_format.h"
#include "culling.h"

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

static void sort_targets(vector<Target>& targets) {
	std::sort(targets.begin(), targets.end(), config.toppercent < 100.0 ? Target::comp_score : Target::comp_evalue);
}

template<typename It>
static It output_range(const It begin, const It end, const Search::Config& cfg) {
	if (end <= begin)
		return begin;
	It i = begin;
	if (i->filter_evalue == DBL_MAX)
		return begin;
	if (config.toppercent < 100.0) {
		const double cutoff = std::max(top_cutoff_score(score_matrix.bitscore(begin->filter_score)), 1.0);
		while (i < end && (score_matrix.bitscore(i->filter_score) >= cutoff))
			++i;
	}
	else {
		i += std::min((ptrdiff_t)cfg.max_target_seqs, end - begin);
		while (--i > begin && i->filter_evalue == DBL_MAX);
		++i;
	}
	return i;
}

bool append_hits(vector<Target>& targets, vector<Target>::iterator begin, vector<Target>::iterator end, bool with_culling, const Search::Config& cfg) {
	if (end <= begin)
		return false;
	bool new_hits = config.toppercent == 100.0 && (int64_t)targets.size() < cfg.max_target_seqs;
	bool append = !with_culling || new_hits;

	culling(targets, append, cfg);
	
	int max_score = 0;
	double min_evalue = DBL_MAX;
	for (auto i = begin; i < end; ++i) {
		max_score = std::max(max_score, i->filter_score);
		min_evalue = std::min(min_evalue, i->filter_evalue);
	}

	vector<Target>::const_iterator range_end = output_range(targets.begin(), targets.end(), cfg);

	if (targets.empty()
		|| (config.toppercent == 100.0 && min_evalue <= (range_end - 1)->filter_evalue)
		|| (config.toppercent != 100.0 && max_score >= top_cutoff_score((range_end - 1)->filter_score))) {
		append = true;
		new_hits = true;
	}

	if(append)
		targets.insert(targets.end(), std::make_move_iterator(begin), std::make_move_iterator(end));

	return new_hits;
}

bool filter_hsp(Hsp& hsp, int source_query_len, const char *query_title, int subject_len, const char* subject_title, const Sequence& query_seq, const Sequence& subject_seq, const double query_self_aln_score, const double target_self_aln_score, const OutputFormat* output_format) {
	bool cluster_threshold = true;
	if (config.cluster_threshold.present()) {
		HspContext context(hsp, 0, 0, TranslatedSequence(query_seq), query_title, 0, subject_len, subject_title, 0, 0, subject_seq, 0, query_self_aln_score, target_self_aln_score);
		RecursiveParser rp(&context, dynamic_cast<const Clustering_format*>(output_format)->format.c_str());
		cluster_threshold = rp.evaluate() >= config.cluster_threshold;
	}
	const double qcov = hsp.query_cover_percent(source_query_len),
		tcov = hsp.subject_cover_percent(subject_len),
		approx_min_id = config.approx_min_id.get(0.0);
	return !cluster_threshold
		|| hsp.id_percent() < config.min_id
		|| (approx_min_id > 0 && hsp.approx_id < approx_min_id)
		|| qcov < config.query_cover
		|| tcov < config.subject_cover
		|| (qcov < config.query_or_target_cover && tcov < config.query_or_target_cover)
		|| (config.no_self_hits
			&& query_seq == subject_seq
			&& strcmp(query_title, subject_title) == 0);
}

void Match::apply_filters(int source_query_len, const char *query_title, const Sequence& query_seq, const double query_self_aln_score, const Block& targets, const OutputFormat* output_format)
{
	const char* title = config.no_self_hits ? targets.ids()[target_block_id] : nullptr;
	const Sequence seq = targets.seqs()[target_block_id];
	const int len = seq.length();
	const double self_aln = targets.has_self_aln() ? targets.self_aln_score(target_block_id) : 0.0;
	for (list<Hsp>::iterator i = hsp.begin(); i != hsp.end();) {
		if (filter_hsp(*i, source_query_len, query_title, len, title, query_seq, seq, query_self_aln_score, self_aln, output_format))
			i = hsp.erase(i);
		else
			++i;
	}
	filter_evalue = hsp.empty() ? DBL_MAX : hsp.front().evalue;
	filter_score = hsp.empty() ? 0 : hsp.front().score;
}

void culling(std::vector<Target>& targets, bool sort_only, const Search::Config& cfg) {
	sort_targets(targets);
	if (!sort_only)
		targets.erase(output_range(targets.begin(), targets.end(), cfg), targets.end());
}

void apply_filters(std::vector<Match>::iterator begin, std::vector<Match>::iterator end, int source_query_len, const char* query_title, const double query_self_aln_score, const Sequence& query_seq, const Search::Config& cfg) {
	if (config.min_id > 0 || config.approx_min_id.get(0.0) > 0 || config.query_cover > 0 || config.subject_cover > 0 || config.query_or_target_cover > 0 || config.no_self_hits || config.cluster_threshold.present())
		for (auto i = begin; i < end; ++i)
			i->apply_filters(source_query_len, query_title, query_seq, query_self_aln_score, *cfg.target, cfg.output_format.get());
}

void culling(std::vector<Match>& targets, const Search::Config& cfg) {
	std::sort(targets.begin(), targets.end(), config.toppercent < 100.0 ? Match::cmp_score : Match::cmp_evalue);
	targets.erase(output_range(targets.begin(), targets.end(), cfg), targets.end());
}

}