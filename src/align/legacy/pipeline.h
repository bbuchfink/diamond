#pragma once
#include "query_mapper.h"

namespace ExtensionPipeline {
	namespace BandedSwipe {
		struct Target;
		struct Pipeline : public QueryMapper
		{
			Pipeline(size_t query_id, Search::Hit* begin, Search::Hit* end, DpStat& dp_stat, const Search::Config& cfg) :
				QueryMapper(query_id, begin, end, cfg),
				dp_stat(dp_stat)
			{}
			Target& target(size_t i);
			virtual void run(Statistics& stat, const Search::Config& cfg) override;
			void run_swipe(bool score_only);
			void range_ranking(const int64_t max_target_seqs);
			DpStat& dp_stat;
		};
	}
}
