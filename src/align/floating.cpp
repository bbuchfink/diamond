// SPDX-License-Identifier: GPL-3.0-or-later
#include "extend.h"
#include "target.h"
#include "culling.h"
#include "data/block/block.h"
#include "dp/scalar/floating_sw.h"
#include "dp/floating_scores.h"

namespace Extension {

std::vector<Match> extend_floating(const Query& query, const Search::Config& cfg) {
	// Intermediate records and DAA store integer raw scores. Reject these paths
	// rather than silently dropping the fractional part of a reference score.
	if (cfg.blocked_processing)
		throw std::runtime_error("--float-scores requires the database to fit in one block; increase --block-size.");
	std::vector<Match> matches;
	const Sequence q = query.sequence.front();
	const auto& bias = query.hauser_correction.front();
	const auto& targets = cfg.target->seqs();
	std::vector<Sequence> target_sequences;
	target_sequences.reserve(targets.size());
	for (BlockId i = 0; i < targets.size(); ++i) target_sequences.push_back(targets[i]);
	std::vector<double> scores;
	if (config.float_score_kernel == "simd")
		scores = DP::floating_scores(q, bias, target_sequences,
			{ score_matrix.matrix32(), score_matrix.gap_open(), score_matrix.gap_extend() });
	else {
		for (Sequence target : target_sequences)
			scores.push_back(DP::floating_sw(q.length(), target.length(),
				score_matrix.gap_open(), score_matrix.gap_extend(),
				[&](size_t i, size_t j) { return double(score_matrix(q[i], target[j])) + double(bias[i]); }).score);
	}
	for (BlockId target_id = 0; target_id < targets.size(); ++target_id) {
		const Sequence target = targets[target_id];
		if (scores[target_id] <= 0.0) continue;
		Hsp hsp(true);
		hsp.floating_score = scores[target_id];
		hsp.bit_score = score_matrix.bitscore(hsp.floating_score, true);
		hsp.evalue = score_matrix.evalue(hsp.floating_score, q.length(), target.length());
		hsp.corrected_bit_score = score_matrix.bitscore_corrected(hsp.floating_score, q.length(), target.length());
		if (config.min_bit_score != 0.0 ? hsp.bit_score < config.min_bit_score : hsp.evalue > config.max_evalue)
			continue;
		Match match(target_id, target, nullptr, 0, 0, hsp.evalue);
		match.hsp.push_back(std::move(hsp));
		matches.push_back(std::move(match));
	}
	// With coordinate/identity filters, a higher-scoring target may later fail.
	// In that case defer top/k culling until all eligible tracebacks are checked.
	const bool alignment_filters = config.min_id > 0 || config.approx_min_id.get(0.0) > 0
		|| config.query_cover > 0 || config.subject_cover > 0 || config.query_or_target_cover > 0
		|| config.no_self_hits || config.cluster_threshold.present();
	if (!alignment_filters) {
		culling(matches, cfg);
		if (cfg.output_format->hsp_values == HspValues::NONE) return matches;
	}
	for (auto it = matches.begin(); it != matches.end();) {
		const Sequence target = it->seq;
		const DP::FloatingAlignment alignment = DP::floating_sw(q.length(), target.length(),
			score_matrix.gap_open(), score_matrix.gap_extend(),
			[&](size_t i, size_t j) { return double(score_matrix(q[i], target[j])) + double(bias[i]); });
		Hsp& hsp = it->hsp.front();
		hsp.query_range = { Loc(alignment.query_begin), Loc(alignment.query_end) };
		hsp.query_source_range = hsp.query_range;
		hsp.subject_range = { Loc(alignment.target_begin), Loc(alignment.target_end) };
		hsp.subject_source_range = hsp.subject_range;
		hsp.target_seq = target;
		size_t i = alignment.query_begin, j = alignment.target_begin;
		char previous = 'M';
		for (char op : alignment.operations) {
			++hsp.length;
			if (op == 'M') {
				if (q[i] == target[j]) {
					hsp.transcript.push_back(op_match);
					++hsp.identities;
				}
				else {
					hsp.transcript.push_back(op_substitution, target[j]);
					++hsp.mismatches;
				}
				if (score_matrix(q[i], target[j]) > 0) ++hsp.positives;
				++i; ++j;
			}
			else {
				if (op != previous) ++hsp.gap_openings;
				++hsp.gaps;
				if (op == 'I') { hsp.transcript.push_back(op_insertion, 1u); ++i; }
				else { hsp.transcript.push_back(op_deletion, target[j]); ++j; }
			}
			previous = op;
		}
		hsp.transcript.push_terminator();
		hsp.approx_id = hsp.id_percent();
		if (filter_hsp(hsp, query.source_length, query.title, target.length(), cfg.target->ids()[it->target_block_id],
			q, target, query.self_alignment_score, 0.0, cfg.output_format.get()))
			it = matches.erase(it);
		else ++it;
	}
	culling(matches, cfg);
	return matches;
}

}
