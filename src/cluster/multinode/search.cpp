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

#include "basic/statistics.h"
#include "multinode.h"
#include "../cluster.h"
#include "run/workflow.h"

using std::vector;
using std::string;
using std::unique_ptr;
using std::tie;
using std::runtime_error;
using std::shared_ptr;

static void run_all_vs_all(Job& job) {
	job.log("Running all-vs-all search for round %d/%d", job.round() + 1, job.round_count());
	const string base_dir = config.tmpdir + PATH_SEPARATOR + "round" + std::to_string(job.round()) + PATH_SEPARATOR;
	config.output_file = base_dir + "alignments.tsv";	
	config.self = true;
	config.query_file.clear();
	config.lin_stage1_query = false;
	config.ext_.clear();
	config.comp_based_stats = 1;
	config.database = job.round() == 0 ? job.root_dir() + "input.faa" :
		job.base_dir(job.round() - 1) + PATH_SEPARATOR + "reps" + PATH_SEPARATOR + "reps_all.faa";
	config.fasta_index_file = job.round() == 0 ? job.root_dir() + "input.faa.faidx" : job.base_dir(job.round() - 1) + PATH_SEPARATOR + "reps" + PATH_SEPARATOR + "reps_all.faa.faidx";
	if (config.db_size == 0)
		throw runtime_error("Database size must be set for cascaded all-vs-all round.");
	tie(config.chunk_size, config.lowmem_) = block_size(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)),
		config.db_size,
		config.sensitivity,
		false,
		config.threads_);
	job.log("Block size: %.2f GB, index chunks: %u", config.chunk_size, config.lowmem_);
	unique_ptr<vector<BitVector>> seed_filter;
	Search::run(seed_filter);
}

static void run_block_combo(Job& job, const VolumedFile& volumes, int64_t r, int64_t i, string base_dir, unique_ptr<vector<BitVector>>& seed_hit_table) {
	config.ext_ = "full";
	config.lin_stage1_query = true;
	if (r == i) {
		config.self = true;
		config.query_file.clear();
	}
	else {
		config.query_file = { volumes[i].path };
		config.self = false;
	}
	config.chunk_size = 1024;
	config.database.clear();
	config.fasta_index_file.clear();
	if (config.db_size == 0)
		throw runtime_error("Database size must be set for cascaded linear search round.");
	config.comp_based_stats = 1; // TODO
	config.output_file = base_dir + std::to_string(r) + "_" + std::to_string(i) + ".tsv";
	log_rss();
	TaskTimer timer("Opening the database");
	shared_ptr<SequenceFile> db(SequenceFile::auto_create({ volumes[r].path }, SequenceFile::Flags::SEQS | SequenceFile::Flags::TITLES | SequenceFile::Flags::RANK_BY_SEQID, amino_acid_traits));
	timer.finish();
	if (!db->open_stats().empty())
		message_stream << db->open_stats();
	Search::run(seed_hit_table, db, nullptr, nullptr);
	job.stats().extensions_computed += statistics.get(Statistics::EXT16) + statistics.get(Statistics::EXT32) + statistics.get(Statistics::EXT8);
}

void run_search(Job& job, const VolumedFile& volumes, int64_t r, int64_t i, string base_dir, unique_ptr<vector<BitVector>>& seed_hit_table) {
	config.command = Config::blastp;
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
	config.db_size = volumes.letter_count();
	config.max_target_seqs_ = 0;
	config.toppercent.unset();
	config.iterate = vector<string>();
	config.iterate.unset();
	config.algo = Config::Algo::DOUBLE_INDEXED;
	config.hamming_dist_boundary_check = true;
	config.mapany = false;
	config.lin_stage1_target = false;
	if(job.is_linear_round())
		run_block_combo(job, volumes, r, i, base_dir, seed_hit_table);
	else
		run_all_vs_all(job);
}