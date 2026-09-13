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
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <cstdarg>
#include <algorithm>
#include <string>
#include "basic/config.h"
#include "volume.h"
#include "multinode.h"
#include "data/sequence_file.h"
#include "cluster/cluster.h"
#include "tools/tools.h"
#include "util/log_stream.h"
#include "util/io/file.h"

const char* const DEFAULT_MEMORY_LIMIT = "16G";
const double CASCADED_ROUND_MAX_EVALUE = 0.001;

using std::runtime_error;
using std::vector;
using std::string;
using std::pair;
using std::unique_ptr;
using std::ofstream;
using std::atomic;
using std::endl;
using std::ifstream;
using std::tie;
using std::shared_ptr;

#ifdef WIN32
static const char* const LOG_EOL = "\r\n";
#else
static const char* const LOG_EOL = "\n";
#endif

std::string Job::log_prefix() const {
	const long long int t = std::chrono::duration_cast<std::chrono::duration<long long int>>(std::chrono::system_clock::now() - start_).count();
	char buffer[128];
	snprintf(buffer, sizeof(buffer), "[%" PRId64 ", %lli] ", worker_id_, t);
	return string(buffer);
}

void Job::log(const char* format, ...) {
	char buffer[1024];
	va_list args;
	va_start(args, format);
	const int i = vsnprintf(buffer, sizeof(buffer), format, args);
	va_end(args);
	if (i < 0)
		return;
	log_raw(string(buffer, std::min((size_t)i, sizeof(buffer) - 1)));
}

void Job::log_raw(const string& message) {
	size_t n = message.size();
	while (n > 0 && (message[n - 1] == '\n' || message[n - 1] == '\r'))
		--n;
	const string s = log_prefix() + message.substr(0, n) + LOG_EOL;
	*message_stream << s;
	log_file_->push(s);
}

void Job::log(const ClusterStats& stats) {
	std::ostringstream ss;
	//stats.masking_stat.print(ss);
	log_raw(ss.str());
	//log("Seeds considered: %" PRIu64, stats.seeds_considered);
	//log("Seeds indexed: %" PRIu64, stats.seeds_indexed);
	log("Extensions computed: %" PRIu64, stats.extensions_computed);
	//log("Alignments passing e-value filter: %" PRIu64, stats.hits_evalue_filtered);
	//log("Alignments passing all filters: %" PRIu64, stats.hits_filtered);
}

static void run_block_combos(Job& job, const VolumedFile& superblocks, const string& base_dir, const string& aln_volumes_path) {
	int64_t r;
	const bool lin_index = use_lin_index(job);
	if (lin_index) {
		configure_round(job, superblocks.letter_count());
		build_lin_indices(job, superblocks);
		if (!job.goon())
			return;
	}
	Atomic q(base_dir + "queue", job);
	Atomic finished(base_dir + "finished", job);
	const int64_t n = (int64_t)superblocks.size();
	while (job.goon() && (r = q.fetch_add(), r < n)) {
		unique_ptr<vector<BitVector>> seed_hit_table;
		if (r > 0)
			seed_hit_table.reset(new vector<BitVector>());
		for (int i = 0; i <= r; ++i) {
			job.log("Searching blocks. Blocks=%lli,%lli", r + 1, i + 1);
			/*if (!seed_hit_table->empty()) {
				for (size_t i = 0; i < seed_hit_table->size(); ++i)
					job.log("Seed hit table paired positions shape %zu: %zu/%zu", i, seed_hit_table->operator[](i).one_count(), seed_hit_table->operator[](i).size());
			}*/
			run_search(job, superblocks, r, i, base_dir, seed_hit_table);
		}
		finished.fetch_add();
		job.finish_step();
	}
	if (!job.goon())
		return;
	finished.await(n);

	Atomic volumes_lock(base_dir + "volumes_lock", job);
	Atomic volumes_done(base_dir + "volumes_done", job);
	if (volumes_lock.fetch_add() == 0) {
		if (lin_index)
			remove_lin_indices(superblocks);
		job.log("Writing alignment volume list");
		ofstream out;
		out.exceptions(std::ios::failbit | std::ios::badbit);
		out.open(aln_volumes_path);
		for (uint64_t r = 0; r < superblocks.size(); ++r) {
			for (uint64_t i = 0; i <= r; ++i) {
				out << base_dir + std::to_string(r) + "_" + std::to_string(i) + ".tsv" << '\n';
			}
		}
		out.close();
		volumes_done.fetch_add();
		job.finish_step();
	}
	else
		volumes_done.await(1);
}

static pair<string, uint64_t> run_round(Job& job, const VolumedFile& superblocks, const string& round_minichunks) {
	if (config.mutual_cover.present()) {
		config.min_length_ratio = config.sensitivity < Sensitivity::LINCLUST_40 ?
			std::min(config.mutual_cover.get_present() / 100 + 0.05, 1.0)
			: config.mutual_cover.get_present() / 100 - 0.05;
	}
	const bool linear = job.is_linear_round();
	job.log("Starting round %i/%i sensitivity=%s linear=%s sequence letters=%" PRIu64, job.round() + 1, job.round_count(), to_string(config.sensitivity).c_str(), linear ? "true" : "false", superblocks.letter_count());
	const int64_t BUF_SIZE = 4096;
	const string base_dir = job.base_dir() + PATH_SEPARATOR + "alignments" + PATH_SEPARATOR;
	const string aln_path = job.base_dir() + "alignments.tsv";
	const string aln_volumes_path = base_dir + "volumes.tsv";
	const bool mutual_cover = config.mutual_cover.present();
	job.make_temp_dir(base_dir);
	if (linear) {
		run_block_combos(job, superblocks, base_dir, aln_volumes_path);
	}
	else {
		unique_ptr<vector<BitVector>> seed_hit_table;
		run_search(job, superblocks, -1, -1, base_dir, seed_hit_table);
	}	
	if (!job.goon())
		return pair<string, uint64_t>("", 0);
	superblocks.remove(false, linear && config.mutual_cover.blank(), false);
	if (job.last_round()) {
		if (!config.fasta_index_file.empty())
			remove_tmp_file(config.fasta_index_file);
		superblocks.remove(job.round() > 0, true, false);
	}	
	Atomic gvc_lock(base_dir + "gvc_lock", job);
	Atomic gvc_done(base_dir + "gvc_done", job);
	GVC::Cfg cfg;
	if (gvc_lock.fetch_add() == 0) {
		job.log("Running greedy vertex cover");
		config.max_oid = job.max_oid();
		config.edges = aln_path;
		config.edge_format = mutual_cover ? "triplet" : "";
		config.symmetric = mutual_cover;
		config.output_file.clear();
		cfg.tmp_dir = job.base_dir();
		if (linear) {
			VolumedFile edges(aln_volumes_path);
			cfg.edges = &edges;
			GVC::greedy_vertex_cover(cfg);
			cfg.edges = nullptr;
			edges.remove(false);
		}
		else {
			GVC::greedy_vertex_cover(cfg);
			remove_tmp_file(aln_path);
		}
		File clusters(job.base_dir() + "clusters.bin", "wb");
		clusters.write(cfg.clustering->data(), cfg.clustering->size() * sizeof(OId));
		clusters.close();
		gvc_done.fetch_add();
		job.finish_step();
	}
	else
		gvc_done.await(1);
	gvc_lock.close();
	gvc_done.close();
	rmdir(base_dir.c_str());
	if (!job.goon())
		return pair<string, uint64_t>("", 0);
	return get_reps(job, round_minichunks, std::move(cfg.clustering));
}

void multinode() {
	if (config.single_step && config.parallel_tmpdir.empty())
		throw runtime_error("Cannot use --single-step without --parallel-tmpdir");
	config.database.require();
	Cluster::init_thresholds();
	const Header hdr_format = TabularFormat::header_format(::Config::cluster);
	const bool parallel = !config.parallel_tmpdir.empty();
	if (config.output_file.empty())
		throw runtime_error("Option missing: output file (--out/-o)");
	const string output_file = config.output_file;
	config.file_buffer_size = 64 * 1024; // TODO
	const bool linclust = config.command == ::Config::LINCLUST;
	const vector<string> rounds = Cluster::cluster_steps(config.approx_min_id.present() ? config.approx_min_id : config.min_id, linclust); // TODO
	if (parallel) {
		for (const string& step : rounds) {
			if (!ends_with(step, "_lin"))
				throw runtime_error("Parallel workflow only supports linclust workflows, support for all-vs-all rounds will be added in a future version.");
		}
	}
	const double evalue_cutoff = config.max_evalue,
		target_approx_id = config.approx_min_id.present() ? config.approx_min_id.get_present() : 0.0;
	const bool anchored_swipe = config.anchored_swipe, is_linclust = Cluster::is_linclust(rounds);
	// TODO
	config.hamming_ext = config.approx_min_id.present() ? config.approx_min_id.get_present() >= 50.0 : false;
	//config.freq_masking = true;
	
	if (parallel) {
		config.tmpdir = config.parallel_tmpdir + PATH_SEPARATOR + "diamond-tmp-" + Const::version_string + PATH_SEPARATOR;
		mkdir(config.tmpdir);
	}
	else
		config.tmpdir = create_temp_directory(config.tmpdir, "diamond-tmp-") + PATH_SEPARATOR;	
	Job job;
	const string input_vols = config.tmpdir + "input_vols.tsv";
	Atomic lock(config.tmpdir + "startup_lock", job), done(config.tmpdir + "startup_done", job);
	if (lock.fetch_add() == 0) {
		ofstream input_vols_file(input_vols);
		input_vols_file << config.database.get_present() << endl;
		input_vols_file.close();
		done.fetch_add();
	}
	else
		done.await(1);	
	config.database = input_vols;
	
	VolumedFile volumes(config.database.get_present());
	if (job.worker_id() == 0) {
		if (config.mutual_cover.present())
			job.log("Bi-directional coverage = %f", config.mutual_cover.get_present());
		else
			job.log("Uni-directional coverage = %f", config.member_cover.get(80));
		if(config.approx_min_id.present())
			job.log("Approximate sequence id cutoff = %f", config.approx_min_id.get(0));
		else
			job.log("Sequence id cutoff = %f", config.min_id);
		job.log("#Volumes = %lli", volumes.size());
	}

	if (max_open_files_per_process() < 1024) {
		const long n = raise_open_files_limit(1024);
		job.log("Raised open files limit to %li", n);
	} else
		job.log("Open files limit = %li", max_open_files_per_process());
	
	string rep_minichunks, input_minichunks_seqs, input_minichunks_accs;
	uint64_t letters = 0;
	job.set_round_count((int)rounds.size(), rounds);
	tie(input_minichunks_seqs, input_minichunks_accs) = len_sort(job, volumes);
	if (!job.goon())
		return;
	const vector<string> ccd_arg = config.connected_component_depth;

	for (size_t i = 0; i < rounds.size(); ++i) {
		const bool linear_round = ends_with(rounds[i], "_lin");
		config.sensitivity = from_string<Sensitivity>(rstrip(rounds[i], "_lin"));
		const vector<string> round_approx_id = config.round_approx_id.empty() ? Cluster::default_round_approx_id(job.round_count()) : config.round_approx_id;
		if (config.min_id == 0.0) {
			config.approx_min_id = std::max(target_approx_id, Cluster::round_value(round_approx_id, "--round-approx-id", i, (int)rounds.size()));
			job.log("Approximate sequence id cutoff (round) = %f", config.approx_min_id.get_present());
		}
		config.max_evalue = i == rounds.size() - 1 ? evalue_cutoff : std::min(evalue_cutoff, CASCADED_ROUND_MAX_EVALUE);
		config.anchored_swipe = anchored_swipe && (linclust || !config.lin_stage1_query);
		if (anchored_swipe)
			config.extension_mode = "banded-fast";
		const int ccd = Cluster::round_ccd(ccd_arg, i, rounds.size(), linear_round);
		config.connected_component_depth.clear();
		config.connected_component_depth.push_back(std::to_string(ccd));

		const string superblocks = job.round() == 0 ? make_merged_blocks(job, input_minichunks_seqs, job.base_dir() + "input_superblocks" + PATH_SEPARATOR, volumes.letter_count())
			: make_merged_blocks(job, rep_minichunks, job.base_dir() + "input_superblocks" + PATH_SEPARATOR, volumes.letter_count());
		VolumedFile input_volumes(superblocks);
		input_volumes.set_letter_count(volumes.letter_count());

		tie(rep_minichunks, letters) = run_round(job, input_volumes, i == 0 ? input_minichunks_seqs : rep_minichunks);
		if (!job.goon())
			return;
		if (i < rounds.size() - 1)
			job.next_round();
#ifdef HAVE_MALLOC_TRIM
		malloc_trim(0);
#endif
	}
	Atomic output_lock(job.root_dir() + PATH_SEPARATOR + "output_lock", job);
	config.output_file = output_file;
	if (output_lock.fetch_add() == 0) {
		VolumedFile acc_vols(input_minichunks_accs);
		merge(job, acc_vols, volumes, hdr_format);
		job.log(job.stats());
		output_lock.close();
		lock.close();
		done.close();
		if (parallel)
			return;
		job.finish();
		remove_tmp_file(input_vols);
		remove_tmp_file(job.root_dir() + "input_minichunks" + PATH_SEPARATOR + "seqs.tsv");
		for (size_t i = 0; i < rounds.size(); ++i) {
			remove_tmp_file(job.base_dir(i) + PATH_SEPARATOR + "reps" + PATH_SEPARATOR + "reps.tsv");
			remove_tmp_file(job.base_dir(i) + PATH_SEPARATOR + "rep_minichunks" + PATH_SEPARATOR + "reps.tsv");
			remove_tmp_file(job.base_dir(i) + PATH_SEPARATOR + "input_superblocks" + PATH_SEPARATOR + "volumes.tsv");
			rmdir(job.base_dir(i) + PATH_SEPARATOR + "input_superblocks");
			rmdir(job.base_dir(i) + PATH_SEPARATOR + "reps");
			rmdir(job.base_dir(i) + PATH_SEPARATOR + "rep_minichunks");
			rmdir(job.base_dir(i));
		}
		//input_volumes.remove(false, false, true);
		rmdir(job.root_dir() + "input_minichunks");
		rmdir(config.tmpdir);
	}
}
