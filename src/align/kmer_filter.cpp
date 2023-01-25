#include "../util/kmer/filter.h"
#include "target.h"
#include "../basic/config.h"
#include "../dp/ungapped.h"

using std::move;
using std::vector;
using std::pair;

namespace Extension {

static const Loc MAX_LEN_DIFF_TRIVIAL_ALN = 3;

pair<SeedHitList, vector<Match>> kmer_filter(Sequence query, const int8_t* query_cbs, const Block& targets, const SeedHitList& l) {
	const KmerFilter filter(query, config.filter_kmer_len);
	SeedHitList r;
	vector<Match> matches;
	for (vector<uint32_t>::const_iterator i = l.target_block_ids.begin(); i != l.target_block_ids.end(); ++i) {
		const Sequence target = targets.seqs()[*i];

		if (abs(query.length() - target.length()) <= MAX_LEN_DIFF_TRIVIAL_ALN) {
			Hsp hsp = trivial(query, target, query_cbs);
			if (hsp.score) {
				matches.emplace_back(*i, target, ::Stats::TargetMatrix(), 0, hsp.score, hsp.evalue);
				matches.back().hsp.push_back(move(hsp));
				matches.back().apply_filters(query.length(), "", query, 0, targets, nullptr);
				continue;
			}
		}

		const pair<double, double> cov = filter.covered(targets.seqs()[*i]);
		if (cov.first >= config.filter_kmer_cutoff || cov.second >= config.filter_kmer_cutoff) {
			const auto target = i - l.target_block_ids.begin();
			r.target_block_ids.push_back(*i);
			r.seed_hits.push_back(l.seed_hits.cbegin(target), l.seed_hits.cend(target));
			r.target_scores.push_back(l.target_scores[target]);
		}
	}
	return { move(r), move(matches) };
}

}