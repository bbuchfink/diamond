#include <memory>
#include "incremental.h"
#include "../basic/config.h"
#include "../../data/sequence_file.h"
#include "../../run/workflow.h"
#include "../../data/block/block_wrapper.h"
#include "../../util/algo/algo.h"
#include "common.h"
#include "../../util/table.h"
#include "../cascaded/cascaded.h"

using std::unique_ptr;
using std::shared_ptr;
using std::vector;
using std::endl;
using std::pair;
using Util::Table;
using Cluster::Callback;

namespace Cluster { namespace Incremental {

using Groups = vector<pair<CentroidId, OId>>;

std::string Algo::get_description() {
	return "Incremental clustering (default)";
}	

struct BestCentroid : public Consumer, public vector<OId> {
	BestCentroid(BlockId block_size):
		vector<OId>(block_size, -1)
	{}
	virtual void consume(const char *ptr, size_t n) {
		const char* end = ptr + n;
		while (ptr < end) {
			const auto edge = *(Output::Format::Edge::Data*)ptr;
			ptr += sizeof(Output::Format::Edge::Data);
			this->operator[](edge.query) = edge.target;
		}
	}
	virtual void finalize() {}
	virtual ~BestCentroid() = default;
};

static void self_align(Block& block, Config& cfg) {
	using Edge = Util::Algo::Edge<SuperBlockId>;
	task_timer timer(("CLUSTER Searching " + std::to_string(block.seqs().size()) + " unaligned sequences").c_str(), cfg.message_stream);
	shared_ptr<Callback> neighbors(new Callback());
	shared_ptr<BlockWrapper> unaligned_wrapper(new BlockWrapper(block));
	config.self = true;
	config.max_target_seqs_ = INT64_MAX;
	config.toppercent = 100;
	config.sensitivity = cfg.sens.back();
	config.chunk_size = 10.0;
	config.mapany = false;
	config.query_or_target_cover = config.member_cover;
	config.query_cover = 0;
	Search::run(unaligned_wrapper, nullptr, neighbors);
	cfg.time_self_aln += timer.seconds();
	const auto n = block.seqs().size();
	cfg.problem_size_self += int64_t(n) * int64_t(n - 1) / 2;

	if (cfg.verbosity >= 2)
		timer.go("CLUSTER Computing clustering");

	message_stream << "Finished search. #Edges: " << neighbors->count << endl;
	timer.go("Allocating buffers");
	vector<Edge> edges(neighbors->count);
	timer.go("Loading edges");
	InputFile f(neighbors->edge_file);
	f.read(edges.data(), neighbors->count);
	f.close_and_delete();
	timer.go("Sorting edges");
	FlatArray<Edge> edge_array = make_flat_array_dense(move(edges), n, config.threads_, Edge::GetKey());
	timer.finish();

	vector<BlockId> c = len_sorted_clust(edge_array); //Util::Algo::greedy_vertex_cover(edge_array, nullptr, false);

	cfg.centroids->init_write();
	int64_t new_centroids = 0;
	vector<CentroidId> block2centroid(n);

	for (BlockId i = 0; i < n; ++i) {
		if (c[i] == i) {
			block2centroid[i] = cfg.centroids->sequence_count();
			cfg.centroids->write_seq(block.seqs()[i], block.ids()[i]);
			cfg.centroid2oid.push_back(block.block_id2oid(i));
			++new_centroids;
		}
	}

	for (BlockId i = 0; i < n; ++i)
		cfg.oid2centroid[block.block_id2oid(i)] = block2centroid[c[i]];

	timer.finish();
	if (cfg.verbosity >= 2)
		cfg.message_stream << "CLUSTER added " << new_centroids << " new centroids, " << cfg.centroids->sequence_count() << " total." << endl;
}

static void search_vs_centroids(Block& block, const int round, Config& cfg) {
	if (cfg.verbosity >= 2)
		cfg.message_stream << "CLUSTER searching vs. centroids sensitivity = " << to_string(cfg.sens[round])
		<< " #sequences = " << block.seqs().size() << " , #centroids = " << cfg.centroids->sequence_count() << endl;
	cfg.status_msg();
	task_timer timer(("Searching " + std::to_string(block.seqs().size()) + " against centroid sequences (" + to_string(cfg.sens[round]) + ")").c_str(), cfg.message_stream);
	shared_ptr<BlockWrapper> block_wrapper(new BlockWrapper(block));
	shared_ptr<BestCentroid> best_centroid(new BestCentroid(block.seqs().size()));
	config.self = false;
	config.max_target_seqs_ = 1;
	config.toppercent = 100;
	config.sensitivity = cfg.sens[round];
	config.chunk_size = std::max(block.seqs().letters() / 1e9 + 0.01, cfg.block_size);
	config.query_or_target_cover = 0;
	config.query_cover = config.member_cover;
	//config.mapany = true;
	cfg.centroids->set_seqinfo_ptr(0);
	Search::run(cfg.centroids, block_wrapper, best_centroid);
	cfg.time_search[round] += timer.seconds();
	cfg.problem_size[round] += block.seqs().size() * cfg.centroids->sequence_count();

	if (cfg.verbosity >= 2)
		timer.go("CLUSTER Assigning to clusters");
	int64_t clustered = 0;
	Block unaligned;
	for (BlockId i = 0; i < block.seqs().size(); ++i) {
		const OId oid = block.block_id2oid(i);
		cfg.oid2centroid[oid] = best_centroid->operator[](i);
		if (best_centroid->operator[](i) == -1)
			unaligned.push_back(block.seqs()[i], block.ids()[i], nullptr, oid, SequenceType::amino_acid, 1);
		else
			++clustered;
	}
	unaligned.seqs().finish_reserve();
	timer.finish();
	if (cfg.verbosity >= 2)
		cfg.message_stream << "CLUSTER " << clustered << " assigned to clusters, " << unaligned.seqs().size() << " unaligned." << endl;
	
	if (round + 1 < cfg.sens.size())
		cfg.cache[round]->append(unaligned);
	else
		self_align(unaligned, cfg);
}

void Algo::run() {
	config.database.require();	
	Config cfg;
	config.db_size = cfg.db->letters();
	if (!config.resume.empty())
		cfg.load_state();
	
	task_timer timer("CLUSTER Opening the input file", cfg.message_stream);
	const int64_t block_size = (int64_t)(cfg.block_size * 1e9), cache_limit = 0; // block_size;
	config.output_format = { "edge" };
	unique_ptr<Block> block;
	if (config.resume.empty()) {
		block.reset(cfg.db->load_seqs(std::min(block_size, config.bootstrap_block)));
		cfg.seqs_processed += block->seqs().size();
		cfg.letters_processed += block->seqs().letters();
		timer.finish();
		self_align(*block, cfg);
	}

	for (;;) {
		timer.go("CLUSTER Loading sequences");
		const int64_t load_size = (int64_t)cfg.centroids->letters() * config.centroid_factor;
		//const int64_t load_size = cfg.letters_processed;
		block.reset(cfg.db->load_seqs(std::min(block_size, load_size)));
		cfg.seqs_processed += block->seqs().size();
		cfg.letters_processed += block->seqs().letters();
		timer.finish();
		if (block->empty())
			break;

		search_vs_centroids(*block, 0, cfg);

		for (int i = 0; i < cfg.cache.size(); ++i)
			if (cfg.cache[i]->seqs().letters() >= std::min(cache_limit, (int64_t)cfg.centroids->letters())) {
				cfg.cache[i]->seqs().finish_reserve();
				search_vs_centroids(*cfg.cache[i], i + 1, cfg);
				cfg.cache[i].reset(new Block);
			}

		if (cfg.verbosity >= 2)
			timer.go("CLUSTER Freeing memory");
		block.reset();

		if (config.timeout && cfg.total_time.seconds() >= config.timeout) {
			cfg.message_stream << "Timeout reached. Next OId = " << cfg.db->tell_seq() << endl;
			cfg.save_state();
			break;
		}
	}

	for (int i = 0; i < cfg.cache.size(); ++i)
		if (cfg.cache[i]->seqs().letters() > 0) {
			cfg.cache[i]->seqs().finish_reserve();
			search_vs_centroids(*cfg.cache[i], i + 1, cfg);
		}

	timer.go("Generating output");
	//const Groups groups = Util::Algo::sort_by_value(cfg.oid2centroid.cbegin(), cfg.oid2centroid.cend(), config.threads_);
	//TextBuffer buf;
	//for (const auto& p : groups) {
//		buf << cfg.centroid2oid[p.first] << '\t' << p.second << '\n';
		//cfg.output_file->write(buf.data(), buf.size());
		//buf.clear();
	//}
	for (int64_t i = 0; i < cfg.oid2centroid.size(); ++i)
		cfg.oid2centroid[i] = cfg.centroid2oid[cfg.oid2centroid[i]];
	output_mem<CentroidId>(*cfg.output_file, *cfg.db, cfg.oid2centroid);

	timer.go("Closing the database");
	cfg.db->close();
	cfg.centroids->close();
	timer.finish();

	Table table;

	table("Total time", cfg.total_time.seconds(), "s");
	table("Self alignment time", cfg.time_self_aln, "s");
	table("Input sequences", cfg.db->sequence_count());
	table("Number of clusters", cfg.centroids->sequence_count());

	for (int i = 0; i < cfg.sens.size(); ++i) {
		table("Time (" + to_string(cfg.sens[i]) + ")", cfg.time_search[i], "s");
		table("Problem size (" + to_string(cfg.sens[i]) + ")", cfg.problem_size[i]);
	}
	table("Problem size self-aln", cfg.problem_size_self);
	cfg.message_stream << endl << table;
}

}}