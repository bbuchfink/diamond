#include "config.h"
#include "../basic/config.h"
#include "../data/block/block.h"
#include "../data/taxonomy_nodes.h"
#include "../data/sequence_file.h"
#include "../search/hit.h"
#include "../util/async_buffer.h"
#include "../search/hit.h"
#include "../util/data_structures/deque.h"
#include "../align/global_ranking/global_ranking.h"
#include "../search/search.h"
#include "../masking/masking.h"
#include "../align/def.h"
#include "../dna/dna_index.h"


using std::endl;
using std::runtime_error;
using std::string;

namespace Search {

Config::Config() :
	self(config.self),
	seed_encoding(config.target_indexed ? SeedEncoding::HASHED : SeedEncoding::SPACED_FACTOR),
	query_masking(MaskingAlgo::NONE),
	target_masking(MaskingAlgo::NONE),
	soft_masking(MaskingAlgo::NONE),
	lazy_masking(false),
	track_aligned_queries(false),
	lin_stage1_target(false),
	max_target_seqs(config.max_target_seqs_.get(25)),
	db(nullptr),
	query_file(nullptr),
	out(nullptr),
	iteration_query_aligned(0)
{
	if (config.iterate.present()) {
		if (config.multiprocessing)
			throw runtime_error("Iterated search is not compatible with --multiprocessing.");
		if (config.target_indexed)
			throw runtime_error("Iterated search is not compatible with --target-indexed.");
		if (config.self)
			throw runtime_error("Iterated search is not compatible with --self.");
		if (config.lin_stage1)
			throw runtime_error("Iterated search is not compatible with --lin-stage1.");
		if (config.linsearch)
			throw runtime_error("Iterated search is not compatible with --linsearch.");
		if (config.iterate.empty()) {
			sensitivity = { {Sensitivity::FASTER, true} };
			const auto rounds = iterated_sens.at(config.sensitivity);
			for (Sensitivity s : rounds)
				sensitivity.push_back(s);
		}
		else
			for (const string& s : config.iterate) {
				if (ends_with(s, "_lin"))
					sensitivity.push_back({ from_string<Sensitivity>(rstrip(s, "_lin")),true });
				else
					sensitivity.push_back(from_string<Sensitivity>(s));
				if (sensitivity.back().sensitivity >= config.sensitivity)
					throw runtime_error("Sensitivity levels set for --iterate must be below target sensitivity.");
			}
	}
	
	sensitivity.push_back(config.sensitivity);
	std::sort(sensitivity.begin(), sensitivity.end());
	if (std::adjacent_find(sensitivity.begin(), sensitivity.end()) != sensitivity.end())
		throw std::runtime_error("The same sensitivity level was specified multiple times for --iterate.");

	if (sensitivity.size() > 1) {
		message_stream << "Running iterated search mode with sensitivity steps:";
		for (Round r : sensitivity)
			message_stream << ' ' << to_string(r.sensitivity);
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

    if(config.command != ::Config::blastn) {
        const MaskingMode masking_mode = from_string<MaskingMode>(config.masking_.get("tantan"));
        switch (masking_mode) {
            case MaskingMode::BLAST_SEG:
                query_masking = MaskingAlgo::NONE;
                target_masking = MaskingAlgo::SEG;
                break;
            case MaskingMode::TANTAN:
                query_masking = MaskingAlgo::TANTAN;
                target_masking = MaskingAlgo::TANTAN;
            default:;
        }
    }
    else {
            if (config.gap_open == -1)
                config.gap_open = 5;
            if (config.gap_extend == -1)
                config.gap_extend = 2;
    }

	if (config.ext_.empty()) {
		if (config.global_ranking_targets || config.swipe_all)
			extension_mode = Extension::Mode::FULL;
		else
			extension_mode = Extension::default_ext_mode.at(sensitivity.back().sensitivity);
	}
	else {
		extension_mode = from_string<Extension::Mode>(config.ext_);
		if (extension_mode != Extension::Mode::FULL) {
			if (config.global_ranking_targets)
				throw std::runtime_error("Global ranking only supports full matrix extension.");
			if (config.swipe_all)
				throw std::runtime_error("--swipe only supports full matrix extension.");
		}
	}

	if (extension_mode == Extension::Mode::FULL) {
		if (config.frame_shift > 0)
			throw std::runtime_error("Frameshift alignment does not support full matrix extension.");
	}

	if (config.freq_masking && config.seed_cut_ != 0.0)
		throw std::runtime_error("Incompatible options: --freq-masking, --seed-cut.");
	if (config.freq_sd_ != 0.0 && !config.freq_masking)
		throw std::runtime_error("--freq-sd requires --freq-masking.");

	if (max_target_seqs == 0)
		max_target_seqs = INT64_MAX;

	if (config.minimizer_window_ && config.algo == ::Config::Algo::CTG_SEED)
		throw runtime_error("Minimizer setting is not compatible with contiguous seed mode.");
}

Config::~Config() {

}

void Config::free()
{
}

}