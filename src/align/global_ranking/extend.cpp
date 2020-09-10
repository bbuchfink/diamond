/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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

#include <unordered_map>
#include <mutex>
#include <memory>
#include <thread>
#include "global_ranking.h"
#include "../output/output.h"
#include "../dp/dp.h"

using std::unique_ptr;
using std::mutex;
using std::thread;
using std::unordered_map;

namespace Extension { namespace GlobalRanking {

typedef unordered_map<uint32_t, uint32_t> TargetMap;

void extend_query(uint32_t query_id, vector<uint32_t>::const_iterator target_begin, vector<uint32_t>::const_iterator target_end, const TargetMap& db2block_id) {
	vector<DpTarget> dp_targets;
	const size_t n = target_end - target_begin;
	dp_targets.reserve(n);
	for (size_t i = 0; i < n; ++i)
		dp_targets.emplace_back(ref_seqs::get()[db2block_id.at(target_begin[i])], i);
}

void align_worker(InputFile* query_list, const TargetMap* db2block_id) {
	pair<uint32_t, vector<uint32_t>> input;
	while (input = fetch_query_targets(*query_list), !input.second.empty()) {
		extend_query(input.first, input.second.begin(), input.second.end(), *db2block_id);
	}
}

void extend(DatabaseFile& db, TempFile& merged_query_list, BitVector& ranking_db_filter, Consumer& master_out) {
	task_timer timer("Loading reference sequences");
	InputFile query_list(merged_query_list);
	vector<uint32_t> block2db_id;
	db.rewind();
 	db.load_seqs(&block2db_id, SIZE_MAX, &ref_seqs::data_, &ref_ids::data_, true, &ranking_db_filter, true);
	TargetMap db2block_id;
	db2block_id.reserve(block2db_id.size());
	for (size_t i = 0; i < block2db_id.size(); ++i)
		db2block_id[block2db_id[i]] = i;
	timer.finish();
	verbose_stream << "#Ranked database sequences: " << ref_seqs::get().get_length() << endl;

	timer.go("Computing alignments");
	OutputSink::instance.reset(new OutputSink(0, &master_out));
	vector<thread> threads;
	for (size_t i = 0; i < config.threads_; ++i)
		threads.emplace_back(align_worker, &query_list, &db2block_id);
	for (auto& i : threads)
		i.join();

	timer.go("Cleaning up");
	query_list.close_and_delete();
	delete ref_seqs::data_;
	delete ref_ids::data_;
}

}}