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
#include <cstdarg>
#include "basic/config.h"
#include "volume.h"
#include "multinode.h"
#include "data/sequence_file.h"
#include "cluster/cluster.h"
#include "tools/tools.h"
#include "util/log_stream.h"

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

void Job::log(const char* format, ...) {
	char buffer[1024];
	const long long int t = std::chrono::duration_cast<std::chrono::duration<long long int>>(std::chrono::system_clock::now() - start_).count();
	char* ptr = buffer + snprintf(buffer, 1024, "[%" PRId64 ", %lli] ", worker_id_, t);
	va_list args;
	va_start(args, format);
	int i = vsnprintf(ptr, 1024 - (ptr - buffer), format, args);
#ifdef WIN32
	ptr[i++] = '\r';
#endif
	ptr[i++] = '\n';
	ptr[i] = '\0';
	*message_stream << buffer;
	log_file_->push(buffer);
	va_end(args);
}

void Job::log(const ClusterStats& stats) {
	std::ostringstream ss;
	//stats.masking_stat.print(ss);
	log(ss.str().c_str());
	//log("Seeds considered: %" PRIu64, stats.seeds_considered);
	//log("Seeds indexed: %" PRIu64, stats.seeds_indexed);
	log("Extensions computed: %" PRIu64, stats.extensions_computed);
	//log("Alignments passing e-value filter: %" PRIu64, stats.hits_evalue_filtered);
	//log("Alignments passing all filters: %" PRIu64, stats.hits_filtered);
}

static void run_block_combos(Job& job, const VolumedFile& volumes, const string& base_dir, const string& aln_path) {
	int64_t r;
	Atomic q(base_dir + "queue", job);
	Atomic finished(base_dir + "finished", job);
	const int64_t n = (int64_t)volumes.size();
	while (r = q.fetch_add(), r < n) {
		unique_ptr<vector<BitVector>> seed_hit_table(new vector<BitVector>());
		for (int i = 0; i <= r; ++i) {
			job.log("Searching blocks. Blocks=%lli,%lli", r + 1, i + 1);
			/*if (!seed_hit_table->empty()) {
				for (size_t i = 0; i < seed_hit_table->size(); ++i)
					job.log("Seed hit table paired positions shape %zu: %zu/%zu", i, seed_hit_table->operator[](i).one_count(), seed_hit_table->operator[](i).size());
			}*/
			run_search(job, volumes, r, i, base_dir, seed_hit_table);
		}
		finished.fetch_add();
	}
	finished.await(n);

	Atomic concat_lock(base_dir + "concat_lock", job);
	Atomic concat_done(base_dir + "concat_done", job);
	if (concat_lock.fetch_add() == 0) {
		job.log("Concatenating alignment files");
		ofstream out(aln_path);
		for (uint64_t r = 0; r < volumes.size(); ++r) {
			for (uint64_t i = 0; i <= r; ++i) {
				const string src = base_dir + std::to_string(r) + "_" + std::to_string(i) + ".tsv";
				std::ifstream in(src, std::ios::binary);
				if (!in.good())
					throw runtime_error("Error opening file " + src);
				if (in.peek() != std::ifstream::traits_type::eof()) {
					out << in.rdbuf();
					if (!out) throw runtime_error("Error writing " + aln_path);
				}
				in.close();
				remove_tmp_file(src);
			}
		}
		out.close();
		concat_done.fetch_add();
	}
	else
		concat_done.await(1);
}

static pair<string, uint64_t> run_round(Job& job, const VolumedFile& volumes) {
	if (config.mutual_cover.present()) {
		config.min_length_ratio = config.sensitivity < Sensitivity::LINCLUST_40 ?
			std::min(config.mutual_cover.get_present() / 100 + 0.05, 1.0)
			: config.mutual_cover.get_present() / 100 - 0.05;
	}
	const bool linear = job.is_linear_round();
	job.log("Starting round %i/%i sensitivity=%s linear=%s", job.round() + 1, job.round_count(), to_string(config.sensitivity).c_str(), linear ? "true" : "false");
	job.set_round(volumes.sparse_records());
	const int64_t BUF_SIZE = 4096;
	const string base_dir = job.base_dir() + PATH_SEPARATOR + "alignments" + PATH_SEPARATOR;
	const string aln_path = job.base_dir() + "alignments.tsv";
	const bool mutual_cover = config.mutual_cover.present();
	job.make_temp_dir(base_dir);
	if(linear)
		run_block_combos(job, volumes, base_dir, aln_path);
	else {
		unique_ptr<vector<BitVector>> seed_hit_table;
		run_search(job, volumes, -1, -1, base_dir, seed_hit_table);
	}
	if (job.last_round()) {
		if (!config.fasta_index_file.empty())
			remove_tmp_file(config.fasta_index_file);
		if (config.reps_out.empty() || job.round() > 0)
			volumes.remove(job.round() > 0);
	}
	Atomic gvc_lock(base_dir + "gvc_lock", job);
	Atomic gvc_done(base_dir + "gvc_done", job);
	if (gvc_lock.fetch_add() == 0) {		
		job.log("Running greedy vertex cover");		
		const string acc_path = Cluster::gvc_input_rep_list(job.round(), job.root_dir(), &job, volumes.max_oid());
		config.database = acc_path;
		config.edges = aln_path;
		config.edge_format = mutual_cover ? "triplet" : "";
		config.symmetric = mutual_cover;
		config.output_file = job.base_dir() + PATH_SEPARATOR + "clusters.tsv";
		GVC::Cfg cfg;
		cfg.tmp_dir = job.base_dir();
		GVC::greedy_vertex_cover(cfg);
		remove_tmp_file(acc_path);
		remove_tmp_file(aln_path);
		if (!job.last_round()) {
			job.log("Writing representative ids");
			ifstream cl(config.output_file);
			OId rep, member;
			size_t vol = 0;
			unique_ptr<ofstream> rep_out(new ofstream(job.base_dir() + "rep_ids" + std::to_string(vol)));
			ofstream rep_out2(job.base_dir() + "rep_ids");
			while (cl >> rep >> member) {
				while (volumes[vol].oid_end <= member) {
					++vol;
					rep_out.reset(new ofstream(job.base_dir() + "rep_ids" + std::to_string(vol)));
				}
				if (rep == member) {
					*rep_out << rep << endl;
					rep_out2 << rep << endl;
				}
			}
			rep_out->close();
		}
		gvc_done.fetch_add();
	}
	else
		gvc_done.await(1);
	rmdir(base_dir.c_str());
	return get_reps(job, volumes);
}

void multinode() {
	//if (config.min_id != 0.0)
		//throw runtime_error("Option not supported for this workflow: --id");
	config.database.require();
	Cluster::init_thresholds();
	const Header hdr_format = TabularFormat::header_format(::Config::cluster);
	const bool parallel = !config.parallel_tmpdir.empty();
	if (config.output_file.empty())
		throw runtime_error("Option missing: output file (--out/-o)");
	const string output_file = config.output_file;
	config.file_buffer_size = 64 * 1024; // TODO
	const bool linclust = config.command == ::Config::LINCLUST;
	const vector<string> steps = Cluster::cluster_steps(config.approx_min_id.present() ? config.approx_min_id : config.min_id, linclust); // TODO
	if (parallel) {
		for (const string& step : steps) {
			if (!ends_with(step, "_lin"))
				throw runtime_error("Parallel workflow only supports linclust workflows, support for all-vs-all rounds will be added in a future version.");
		}
	}
	const double evalue_cutoff = config.max_evalue,
		target_approx_id = config.approx_min_id.present() ? config.approx_min_id.get_present() : 0.0;
	const bool anchored_swipe = config.anchored_swipe, is_linclust = Cluster::is_linclust(steps);
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

	if (max_open_files_per_process() < 1024)
		raise_open_files_limit(1024);
	
	string reps;
	uint64_t letters = 0;
	job.set_round_count((int)steps.size(), steps);
	const string input_parts = len_sort(job, volumes);
	VolumedFile input_volumes(input_parts);
	input_volumes.set_letter_count(volumes.letter_count());
	const vector<string> ccd_arg = config.connected_component_depth;

	for (size_t i = 0; i < steps.size(); ++i) {
		const bool linear_round = ends_with(steps[i], "_lin");
		config.sensitivity = from_string<Sensitivity>(rstrip(steps[i], "_lin"));
		const vector<string> round_approx_id = config.round_approx_id.empty() ? Cluster::default_round_approx_id(job.round_count()) : config.round_approx_id;
		if (config.min_id == 0.0) {
			config.approx_min_id = std::max(target_approx_id, Cluster::round_value(round_approx_id, "--round-approx-id", i, (int)steps.size()));
			job.log("Approximate sequence id cutoff (round) = %f", config.approx_min_id.get_present());
		}
		config.max_evalue = i == steps.size() - 1 ? evalue_cutoff : std::min(evalue_cutoff, CASCADED_ROUND_MAX_EVALUE);
		config.anchored_swipe = anchored_swipe && (linclust || !config.lin_stage1_query);
		if (anchored_swipe)
			config.ext_ = "banded-fast";
		const int ccd = Cluster::round_ccd(ccd_arg, i, steps.size(), linear_round);
		config.connected_component_depth.clear();
		config.connected_component_depth.push_back(std::to_string(ccd));
		tie(reps, letters) = run_round(job, i == 0 ? input_volumes : VolumedFile(reps, letters));
		if (i < steps.size() - 1)
			job.next_round();
	}
	Atomic output_lock(job.root_dir() + PATH_SEPARATOR + "output_lock", job);
	config.output_file = output_file;
	if (output_lock.fetch_add() == 0) {
		const vector<OId> merged = build_merged(job);
		if (!config.reps_out.empty())
			write_representatives(job, input_volumes, merged);
		merge(job, input_volumes, hdr_format, merged);
		job.log(job.stats());
		output_lock.close();
		lock.close();
		done.close();
		job.finish();
		remove_tmp_file(input_vols);
		rmdir(config.tmpdir);
	}
}