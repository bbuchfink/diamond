/****
DIAMOND protein aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <inttypes.h>
#include <algorithm>
#include "basic/statistics.h"
#include "masking/def.h"
#include "multinode.h"
#include "cluster.h"
#include "run/workflow.h"
#include "util/log_stream.h"
#include "util/system/system.h"

using std::vector;
using std::string;
using std::unique_ptr;
using std::tie;
using std::runtime_error;
using std::shared_ptr;

static const int LINCLUST_QUERY_BINS = 128;

/* Masking algorithm applied to the input when the length sorted minichunks are written.
   Resolved once from the user setting, as configure_round() clears config.masking_ for
   the alignment workflow. Clustering aligns the input against itself, so the asymmetric
   BLAST_SEG mode (target only) is equivalent to masking everything. */
MaskingAlgo input_masking_algo() {
	static const MaskingAlgo algo = []() {
		const MaskingMode mode = config.masking_.present() ? from_string<MaskingMode>(config.masking_.get_present()) : MaskingMode::BLAST_SEG_ALL;
		switch (mode) {
		case MaskingMode::NONE:
			return MaskingAlgo::NONE;
		case MaskingMode::TANTAN:
			return MaskingAlgo::TANTAN;
		default:
			return MaskingAlgo::SEG;
		}
	}();
	return algo;
}

static void run_all_vs_all(Job& job) {
	job.log("Running all-vs-all search for round %d/%d", job.round() + 1, job.round_count());
	const string base_dir = config.tmpdir + PATH_SEPARATOR + "round" + std::to_string(job.round()) + PATH_SEPARATOR;
	config.output_file = base_dir + "alignments.tsv";	
	config.self = true;
	config.query_file.clear();
	config.query_bins.unset();
	config.lin_stage1_query = false;
	config.gapped_filter_evalue_ = -1.0;
	config.anchored_swipe = false;
	//config.comp_based_stats = 1;
	config.database = job.round() == 0 ? job.root_dir() + "input_minichunks" + PATH_SEPARATOR + "0.faa" :
		job.base_dir(job.round() - 1) + PATH_SEPARATOR + "rep_minichunks" + PATH_SEPARATOR + "reps_all.faa";
	config.fasta_index_file = job.round() == 0 ? job.root_dir() + "input.faa.faidx" : job.base_dir(job.round() - 1) + PATH_SEPARATOR + "rep_minichunks" + PATH_SEPARATOR + "reps_all.faa.faidx";
	if (config.db_size == 0)
		throw runtime_error("Database size must be set for cascaded all-vs-all round.");
	tie(config.chunk_size, config.lowmem_) = block_size(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)),
		config.db_size,
		config.sensitivity,
		false,
		config.threads_, config.mutual_cover.present());
	job.log("Block size: %.2f GB, index chunks: %u", config.chunk_size, config.lowmem_);
	unique_ptr<vector<BitVector>> seed_filter;
	Search::run(seed_filter);
}

static void run_block_combo(Job& job, const VolumedFile& volumes, int64_t r, int64_t i, string base_dir, unique_ptr<vector<BitVector>>& seed_hit_table) {
	config.lin_stage1_query = true;	
	config.lin_index_file = use_lin_index(job) ? lin_index_file(volumes[i].path) : string();
	if (r == i) {
		config.self = true;
		config.query_file.clear();
	}
	else {
		config.query_file = { volumes[i].path };
		config.self = false;
	}
	config.gapped_filter_evalue_ = 0.0;
	config.chunk_size = 65536;
	config.query_bins.set_if_blank(LINCLUST_QUERY_BINS);
	config.anchored_swipe = config.comp_based_stats_.get_present() == 0;
	config.database.clear();
	config.fasta_index_file.clear();
	if (config.db_size == 0)
		throw runtime_error("Database size must be set for cascaded linear search round.");
	//config.comp_based_stats = 1; // TODO
	config.output_file = base_dir + std::to_string(r) + "_" + std::to_string(i) + ".tsv";
	log_rss();
	TaskTimer timer("Opening the database");
	shared_ptr<SequenceFile> db;
	try {
		db.reset(SequenceFile::auto_create({ volumes[r].path }, SequenceFile::Flags::SEQS | SequenceFile::Flags::TITLES | SequenceFile::Flags::RANK_BY_SEQID, amino_acid_traits));
	}
	catch (FormatDetectionError& e) {
		throw runtime_error(string("Error opening database chunk file " + volumes[r].path + ": ") + e.what());
	}
	timer.finish();
	if (!db->open_stats().empty())
		*message_stream << db->open_stats();
	Search::run(seed_hit_table, db, nullptr, nullptr);
	job.stats().extensions_computed += statistics.get(Statistics::EXT16) + statistics.get(Statistics::EXT32) + statistics.get(Statistics::EXT8);
}

void configure_round(Job& job, uint64_t letter_count) {
	static const double MATRIX_ADJUST_MAX_ID = 50;
	config.command = Config::blastp;
	config.lin_index_file.clear();
	const bool mutual_cover = config.mutual_cover.present();
	const vector<string> round_coverage = config.round_coverage.empty() ? Cluster::default_round_cov(job.round_count()) : config.round_coverage;
	const double cov_cutoff = mutual_cover ? config.mutual_cover.get_present() : config.member_cover,
		round_cov_cutoff = std::max(cov_cutoff, Cluster::round_value(round_coverage, "--round-coverage", job.round(), job.round_count()));
	if (mutual_cover) {
		config.query_cover = config.subject_cover = round_cov_cutoff;
	}
	else {
		config.query_cover = 0;
		config.subject_cover = 0;
		config.query_or_target_cover = round_cov_cutoff;
	}
	config.output_format = mutual_cover ? vector<string> { "tab", "qseqid", "sseqid", "corrected_bitscore" } : vector<string>{ "tab", "qseqid", "sseqid", "qcovhsp", "scovhsp", "corrected_bitscore" };
	statistics.reset();
	if (letter_count > 0)
		config.db_size = letter_count;
	job.log("Database letter count: %" PRIu64 " maximum OId: %" PRIu64, config.db_size, job.max_oid());
	config.max_target_seqs_ = 0;
	config.toppercent.unset();
	config.iterate = vector<string>();
	if (config.comp_based_stats_.blank()) {
		if (config.approx_min_id.present() && config.approx_min_id.get_present() >= MATRIX_ADJUST_MAX_ID || config.min_id >= MATRIX_ADJUST_MAX_ID)
			config.comp_based_stats_ = 0;
		else
			config.comp_based_stats_ = 6;
	}
	// The input is hard masked when the length sorted minichunks are written
	// (see input_masking_algo), so the alignment workflow does not mask again.
	config.masking_ = "none";
	config.iterate.unset();
	config.algo = Config::Algo::DOUBLE_INDEXED;
	config.mapany = false;
	config.lin_stage1_target = false;
	config.symmetrize_evalue = true;
	config.no_reorder = true;
	config.no_mem_pool = true;
	config.oid_title_max = job.max_oid();
	//config.ungapped_filter_query_len = 1;
	config.output_header.clear();
	config.output_header.unset();
}

void run_search(Job& job, const VolumedFile& volumes, int64_t r, int64_t i, string base_dir, unique_ptr<vector<BitVector>>& seed_hit_table) {
	configure_round(job, volumes.letter_count());
	if(job.is_linear_round())
		run_block_combo(job, volumes, r, i, base_dir, seed_hit_table);
	else
		run_all_vs_all(job);
}