#include "config.h"
#include "../basic/config.h"
#include "../data/block.h"
#include "../data/taxonomy_nodes.h"
#include "../data/sequence_file.h"
#include "../search/hit.h"
#include "../util/async_buffer.h"
#include "../search/hit.h"
#include "../util/data_structures/deque.h"
#include "../align/global_ranking/global_ranking.h"
#include "../search/search.h"

using std::endl;

namespace Search {

Config::Config() :
	self(config.self),
	seed_encoding(config.target_indexed ? SeedEncoding::HASHED : SeedEncoding::SPACED_FACTOR),
	lazy_masking(false),
	track_aligned_queries(false),
	db(nullptr),
	query_file(nullptr),
	out(nullptr),
	taxon_nodes(nullptr),
	taxonomy_scientific_names(nullptr),
	iteration_query_aligned(0)
{
	if (!config.iterate.empty()) {
		if (config.multiprocessing)
			throw std::runtime_error("Iterated search is not compatible with --multiprocessing.");
		if (config.target_indexed)
			throw std::runtime_error("Iterated search is not compatible with --target-indexed.");
		//sensitivity = iterated_sens.find(config.sensitivity)->second;
		for (const string& s : config.iterate)
			sensitivity.push_back(from_string<Sensitivity>(s));
	}
	
	sensitivity.push_back(config.sensitivity);

	if (sensitivity.size() > 1) {
		message_stream << "Running iterated search mode with sensitivity steps:";
		for (Sensitivity s : sensitivity)
			message_stream << ' ' << to_string(s);
		message_stream << endl;
		track_aligned_queries = true;
	}

	if (!config.unaligned.empty() || !config.aligned_file.empty())
		track_aligned_queries = true;

	if (config.multiprocessing && (!config.taxonlist.empty() || !config.taxon_exclude.empty()))
		throw std::runtime_error("Multiprocessing mode is not compatible with database filtering.");

	if (config.global_ranking_targets) {
		if (config.frame_shift)
			throw std::runtime_error("Global ranking mode is not compatible with frameshift alignments.");
		if(config.multiprocessing)
			throw std::runtime_error("Global ranking mode is not compatible with --multiprocessing.");
	}

	if (config.target_indexed && config.algo != ::Config::Algo::AUTO && config.algo != ::Config::Algo::DOUBLE_INDEXED)
		throw std::runtime_error("--target-indexed requires --algo 0");
}

Config::~Config() {

}

void Config::free()
{
	delete taxon_nodes;
	delete taxonomy_scientific_names;
	taxon_nodes = nullptr;
	taxonomy_scientific_names = nullptr;
}

}