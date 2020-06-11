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
#include "../search/trace_pt_buffer.h"
#include "../util/task_queue.h"
#include "../basic/statistics.h"
#include "legacy/query_mapper.h"
#include "../data/metadata.h"

struct Output_writer
{
	Output_writer(OutputFile* f) :
		f_(f)
	{ }
	void operator()(TextBuffer &buf)
	{
		f_->write(buf.get_begin(), buf.size());
		buf.clear();
	}
private:
	OutputFile* const f_;
};

void align_queries(Trace_pt_buffer &trace_pts, Consumer* output_file, const Parameters &params, const Metadata &metadata);

namespace ExtensionPipeline {
	namespace Swipe {
		struct Pipeline : public QueryMapper
		{
			Pipeline(const Parameters &params, size_t query_id, hit* begin, hit* end, const Metadata &metadata) :
				QueryMapper(params, query_id, begin, end, metadata)
			{}
			virtual void run(Statistics &stat, const sequence *subjects = nullptr, size_t subject_count = 0);
			virtual ~Pipeline() {}
		};
	}
	namespace BandedSwipe {
		struct Target;
		struct Pipeline : public QueryMapper
		{
			Pipeline(const Parameters &params, size_t query_id, hit* begin, hit* end, DpStat &dp_stat, const Metadata &metadata, bool target_parallel) :
				QueryMapper(params, query_id, begin, end, metadata, target_parallel),
				dp_stat(dp_stat)
			{}
			Target& target(size_t i);
			virtual void run(Statistics &stat, const sequence *subjects = nullptr, size_t subject_count = 0);
			void run_swipe(bool score_only);
			void range_ranking();
			DpStat &dp_stat;
		};
	}
}