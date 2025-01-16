/****
DIAMOND protein aligner
Copyright (C) 2016-2023 Max Planck Society for the Advancement of Science e.V.
						Benjamin Buchfink

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#include <numeric>
#include "cascaded.h"
#include "util/algo/algo.h"
#include "basic/statistics.h"
#include "run/workflow.h"
#include "util/log_stream.h"

const char* const DEFAULT_MEMORY_LIMIT = "16G";
const double CASCADED_ROUND_MAX_EVALUE = 0.001;

using std::vector;
using std::shared_ptr;
using std::endl;
using std::runtime_error;
using std::numeric_limits;
using std::to_string;
using std::tie;
using std::string;
using std::pair;
using std::iota;
using std::make_pair;

namespace Cluster {

string Cascaded::get_description() {
	return "Cascaded greedy vertex cover algorithm";
}

static BitVector rep_bitset(const vector<SuperBlockId> &centroid, const BitVector *superset = nullptr) {
	BitVector r(centroid.size());
	for (SuperBlockId c : centroid)
		if (!superset || superset->get(c))
			r.set(c);
	return r;
}

vector<SuperBlockId> cluster(shared_ptr<SequenceFile>& db, const shared_ptr<BitVector>& filter, const SuperBlockId* member_counts, int round, int round_count) {
	using Edge = Util::Algo::Edge<SuperBlockId>;
	statistics.reset();
	const bool mutual_cover = config.mutual_cover.present();
	config.command = Config::blastp;
	config.output_format = { "edge" };
	const vector<string> round_coverage = config.round_coverage.empty() ? default_round_cov(round_count) : config.round_coverage;
	const double cov_cutoff = config.mutual_cover.present() ? config.mutual_cover.get_present() : config.member_cover,
		round_cov_cutoff = std::max(cov_cutoff, round_value(round_coverage, "--round-coverage", round, round_count));
	if (config.mutual_cover.present()) {
		config.query_cover = config.subject_cover = round_cov_cutoff;
	}
	else {
		config.query_cover = 0;
		config.subject_cover = 0;
		config.query_or_target_cover = round_cov_cutoff;
	}
	config.algo = Config::Algo::DOUBLE_INDEXED;
	config.max_target_seqs_ = INT64_MAX;
	config.self = true;
	config.iterate.unset();
	config.mapany = false;
	config.linsearch = false;
	//if (config.db_size == 0 && filter) {
	if(filter) {
		config.db_size = db->letters_filtered(*filter);
	}
	tie(config.chunk_size, config.lowmem_) = config.lin_stage1 && round == 0 ? make_pair<double, int>(32768, 1) :
		block_size(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)),
			filter ? config.db_size : db->letters(),
			config.sensitivity,
			config.lin_stage1,
			config.threads_);

	shared_ptr<Callback> callback(mutual_cover ? (Callback*)new CallbackBidirectional : (Callback*)new CallbackUnidirectional);

	Search::run(db, nullptr, callback, filter);

	message_stream << "Finished search. #Edges: " << callback->count << endl;
	TaskTimer timer("Allocating buffers");
	vector<Edge> edges(callback->count);
	timer.go("Loading edges");
	InputFile f(callback->edge_file);
	f.read(edges.data(), callback->count);
	f.close_and_delete();
	if (!config.aln_out.empty())
		output_edges(config.aln_out, *db, edges);
	timer.go("Sorting edges");
	db->reopen();
	FlatArray<Edge> edge_array = make_flat_array_dense(std::move(edges), (SuperBlockId)db->sequence_count(), config.threads_, Edge::GetKey());
	timer.finish();

	const auto algo = from_string<GraphAlgo>(config.graph_algo);
	const bool last_round = round == round_count - 1;
	const int ccd = round_ccd(round, round_count);
	return algo == GraphAlgo::GREEDY_VERTEX_COVER ?
		Util::Algo::greedy_vertex_cover(edge_array, config.weighted_gvc ? member_counts : nullptr, last_round && !config.strict_gvc && !mutual_cover, !config.no_gvc_reassign, ccd)
		: len_sorted_clust(edge_array);
}

static pair<vector<SuperBlockId>, BitVector> update_clustering(const BitVector& previous_filter, const vector<SuperBlockId>& previous_centroids, vector<SuperBlockId>&& current_centroids, int round) {
	BitVector oid_filter(round == 0 ? rep_bitset(current_centroids) : rep_bitset(current_centroids, &previous_filter));
	if (round > 0)
		for (size_t i = 0; i < current_centroids.size(); ++i)
			if (!previous_filter.get(i))
				current_centroids[i] = current_centroids[previous_centroids[i]];
	return { current_centroids, oid_filter };
}

vector<SuperBlockId> cascaded(shared_ptr<SequenceFile>& db, bool linear) {
	if (db->sequence_count() > (int64_t)numeric_limits<SuperBlockId>::max())
		throw runtime_error("Workflow supports a maximum of " + to_string(numeric_limits<SuperBlockId>::max()) + " input sequences.");
	const auto steps = cluster_steps(config.approx_min_id, linear);
	const double evalue_cutoff = config.max_evalue,
		target_approx_id = config.approx_min_id;
	shared_ptr<BitVector> oid_filter(new BitVector);
	int64_t cluster_count = db->sequence_count();
	vector<SuperBlockId> centroids(cluster_count);
	iota(centroids.begin(), centroids.end(), 0);

	for (int i = 0; i < (int)steps.size(); i++) {
		TaskTimer timer;
		config.lin_stage1 = ends_with(steps[i], "_lin");
		config.sensitivity = from_string<Sensitivity>(rstrip(steps[i], "_lin"));
		if (config.lin_stage1 && config.sensitivity > Sensitivity::SHAPES30x10)
			throw std::runtime_error("Linearization is only supported for the faster, fast, and shapes-30x10 sensitivity modes.");
		const vector<string> round_approx_id = config.round_approx_id.empty() ? default_round_approx_id((int)steps.size()) : config.round_approx_id;
		config.approx_min_id = std::max(target_approx_id, round_value(round_approx_id, "--round-approx-id", i, (int)steps.size()));
		config.max_evalue = (size_t)i == steps.size() - 1 ? evalue_cutoff : std::min(evalue_cutoff, CASCADED_ROUND_MAX_EVALUE);
		tie(centroids, *oid_filter) = update_clustering(*oid_filter,
			centroids,
			cluster(db, i == 0 ? nullptr : oid_filter, config.weighted_gvc ? member_counts(centroids).data() : nullptr, i, (int)steps.size()),
			i);
		const int64_t n = oid_filter->one_count();
		message_stream << "Clustering round " << i + 1 << " complete. #Input sequences: " << cluster_count
			<< " #Clusters: " << n
			<< " #Letters: " << db->letters_filtered(*oid_filter)
			<< " Time: " << timer.seconds() << 's' << endl;
		cluster_count = n;
	}
	return centroids;
}
	
}