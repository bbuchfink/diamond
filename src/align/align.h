/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include <memory>
#include <vector>
#include <map>
#include "../basic/statistics.h"
#include "legacy/query_mapper.h"
#include "../run/workflow.h"
#include "../search/hit.h"

struct Output_writer
{
	Output_writer(OutputFile* f) :
		f_(f)
	{ }
	void operator()(TextBuffer &buf)
	{
		f_->write(buf.data(), buf.size());
		buf.clear();
	}
private:
	OutputFile* const f_;
};

void align_queries(Consumer* output_file, Search::Config &cfg);

namespace ExtensionPipeline {
	namespace Swipe {
		struct Pipeline : public QueryMapper
		{
			Pipeline(size_t query_id, Search::Hit* begin, Search::Hit* end, const Search::Config &cfg) :
				QueryMapper(query_id, begin, end, cfg)
			{}
			virtual void run(Statistics &stat, const Search::Config& cfg) override;
			virtual ~Pipeline() {}
		};
	}
	namespace BandedSwipe {
		struct Target;
		struct Pipeline : public QueryMapper
		{
			Pipeline(size_t query_id, Search::Hit* begin, Search::Hit* end, DpStat &dp_stat, const Search::Config &cfg) :
				QueryMapper(query_id, begin, end, cfg),
				dp_stat(dp_stat)
			{}
			Target& target(size_t i);
			virtual void run(Statistics &stat, const Search::Config& cfg) override;
			void run_swipe(bool score_only);
			void range_ranking(const int64_t max_target_seqs);
			DpStat &dp_stat;
		};
	}
}
