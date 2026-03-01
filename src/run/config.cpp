/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#include "config.h"
#include "basic/config.h"
#include "data/block/block.h"
#include "data/taxonomy_nodes.h"
#include "data/sequence_file.h"
#include "search/hit.h"
#include "search/hit_buffer.h"
#include "search/hit.h"
#include "util/data_structures/deque.h"
#include "align/global_ranking/global_ranking.h"
#include "search/search.h"
#include "masking/masking.h"
#include "align/def.h"

#ifdef WITH_DNA
#include "../dna/extension.h"
#include "../dna/timer.h"
#include "../dna/dna_index.h"
#endif

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
	max_target_seqs(0),
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
		if (config.lin_stage1_query)
			throw runtime_error("Iterated search is not compatible with --lin-stage1.");
		if (config.iterate.empty()) {
			sensitivity = { {Sensitivity::FASTER, true} };
			const auto rounds = iterated_sens.at(config.sensitivity);
			sensitivity.insert(sensitivity.end(), rounds.begin(), rounds.end());
		}
		else {
			const Round target = Round(config.sensitivity, config.lin_stage1_query || config.lin_stage1_target);
			for (const string& s : config.iterate) {
				if (ends_with(s, "_lin"))
					sensitivity.push_back({ from_string<Sensitivity>(rstrip(s, "_lin")),true });
				else
					sensitivity.push_back(from_string<Sensitivity>(s));
				if (!(sensitivity.back() < target))
					throw runtime_error("Sensitivity levels set for --iterate must be below target sensitivity.");
			}
		}
	}
	
	if (sensitivity.empty() || (!sensitivity.empty() && sensitivity.back() != Round(config.sensitivity, config.lin_stage1_target)))
		sensitivity.emplace_back(config.sensitivity, config.lin_stage1_target);
	std::sort(sensitivity.begin(), sensitivity.end());
	if (std::adjacent_find(sensitivity.begin(), sensitivity.end()) != sensitivity.end())
		throw std::runtime_error("The same sensitivity level was specified multiple times for --iterate.");

	if (sensitivity.size() > 1) {
		message_stream << "Running iterated search mode with sensitivity steps: ";
		for (auto it = sensitivity.begin(); it != sensitivity.end(); ++it) {
			message_stream << to_string(it->sensitivity);
			if (it->linearize)
				message_stream << " (linear)";
			if (it != sensitivity.end() - 1)
				message_stream << ", ";
		}
		message_stream << endl;
		track_aligned_queries = true;
	}

	if (!config.unaligned.empty() || !config.aligned_file.empty())
		track_aligned_queries = true;

	if (config.multiprocessing && (!config.taxonlist.empty() || !config.taxon_exclude.empty()))
		throw runtime_error("Multiprocessing mode is not compatible with database filtering.");

	if (config.global_ranking_targets) {
		if (config.frame_shift)
			throw runtime_error("Global ranking mode is not compatible with frameshift alignments.");
		if(config.multiprocessing)
			throw runtime_error("Global ranking mode is not compatible with --multiprocessing.");
	}

	if (config.target_indexed && config.algo != ::Config::Algo::AUTO && config.algo != ::Config::Algo::DOUBLE_INDEXED)
		throw runtime_error("--target-indexed requires --algo 0");

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
#ifdef WITH_DNA
            timer.reset(new Dna::TotalTime());
#endif

    }

	if (config.freq_masking && config.seed_cut_ != 0.0)
		throw runtime_error("Incompatible options: --freq-masking, --seed-cut.");
	if (config.freq_sd_ != 0.0 && !config.freq_masking)
		throw runtime_error("--freq-sd requires --freq-masking.");

	if (config.minimizer_window_ && config.algo == ::Config::Algo::CTG_SEED)
		throw runtime_error("Minimizer setting is not compatible with contiguous seed mode.");

	if (config.query_cover >= 50 && config.query_cover == config.subject_cover && config.min_length_ratio == 0.0 && !align_mode.query_translated) {
		min_length_ratio = config.lin_stage1_query && sensitivity.back().sensitivity < Sensitivity::LINCLUST_40
			? std::min(config.query_cover / 100 + 0.05, 0.92) : std::max(config.query_cover / 100 - 0.05, 0.0);
	}
	else {
		if (align_mode.query_translated && config.min_length_ratio != 0.0)
			throw runtime_error("--min-len-ratio is not supported for translated searches");
		min_length_ratio = config.min_length_ratio;
	}
	log_stream << "Min length ratio: " << min_length_ratio << endl;
	output_format.reset(init_output(max_target_seqs));
}

Config::~Config() {

}

void Config::free()
{
}

}