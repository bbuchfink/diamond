#include "../search/search.h"
#include "load_hits.h"

namespace Extension {

static const int64_t CAP = 1000;

TextBuffer* pipeline_short(BlockId query, Search::Hit* begin, Search::Hit* end, Search::Config& cfg, Statistics& stats) {
	TextBuffer* out = nullptr;
	SeedHitList l = load_hits(begin, end, cfg.target->seqs());
	stats.inc(Statistics::TARGET_HITS0, l.target_block_ids.size());
	std::sort(l.target_scores.begin(), l.target_scores.end());
	stats.inc(Statistics::TARGET_HITS1, std::min((int64_t)l.target_scores.size(), CAP));
	return out;
}

}