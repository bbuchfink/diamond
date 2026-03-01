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

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
#include <inttypes.h>
#include <cstdarg>
#include "basic/config.h"
#include "volume.h"
#include "cluster/cascaded/cascaded.h"
#include "multinode.h"
#include "run/workflow.h"
#include "data/sequence_file.h"
#include "util/sequence/sequence.h"
#include "basic/statistics.h"
#include "cluster/cascaded/cascaded.h"
#include "cluster/cluster.h"
#include "data/fasta/fasta_file.h"

void greedy_vertex_cover();

using std::runtime_error;
using std::vector;
using std::string;
using std::pair;
using std::unique_ptr;
using std::ofstream;
using std::atomic;
using std::endl;
using std::ifstream;

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
	log_stream << buffer;
	log_file_->push(buffer);
	va_end(args);
}

void Job::log(const ClusterStats& stats) {
	std::ostringstream ss;
	stats.masking_stat.print(ss);
	log(ss.str().c_str());
	log("Seeds considered: %" PRIu64, stats.seeds_considered);
	log("Seeds indexed: %" PRIu64, stats.seeds_indexed);
	log("Extensions computed: %" PRIu64, stats.extensions_computed);
	log("Alignments passing e-value filter: %" PRIu64, stats.hits_evalue_filtered);
	log("Alignments passing all filters: %" PRIu64, stats.hits_filtered);
}

static int64_t combos(int64_t n) {
	return n * (n + 1) / 2;
}

static int64_t combo_to_rank(int64_t i, int64_t j, int64_t n) {
	return (2 * n - i + 1) * i / 2 + j - i;
}

static pair<int64_t, int64_t> rank_to_combo(int64_t r, int64_t n) {
	const int64_t a = 2 * n + 1,
		i = std::floor((2 * n + 1 - sqrt(a * a - 8 * r)) / 2),
		j = i + r - (2 * n - i + 1) * i / 2;
	return { i,j };
}

static void run_block_combo(Job& job, const VolumedFile& volumes, int64_t bi, int64_t bj, string base_dir) {
	statistics.reset();
	const bool mutual_cover = config.mutual_cover.present();
	config.command = Config::blastp;
	const vector<string> round_coverage = config.round_coverage.empty() ? Cluster::default_round_cov(job.round_count()) : config.round_coverage;
	const double cov_cutoff = config.mutual_cover.present() ? config.mutual_cover.get_present() : config.member_cover,
		round_cov_cutoff = std::max(cov_cutoff, Cluster::round_value(round_coverage, "--round-coverage", job.round(), job.round_count()));
	if (config.mutual_cover.present()) {
		config.query_cover = config.subject_cover = round_cov_cutoff;
	}
	else {
		config.query_cover = 0;
		config.subject_cover = 0;
		config.query_or_target_cover = round_cov_cutoff;
	}	
	config.self = false;
	config.max_target_seqs_ = 0;
	config.toppercent.unset();
	if (config.mutual_cover.present()) {
		config.query_cover = config.subject_cover = config.mutual_cover.get_present();
	}
	else {
		config.query_cover = config.member_cover;
		config.subject_cover = 0;
	}
	config.query_or_target_cover = 0;
	config.iterate = vector<string>();
	config.iterate.unset();
	if (job.round() == 0) {
		config.qnum_offset = volumes[bi].oid_begin;
		config.snum_offset = volumes[bj].oid_begin;
		config.output_format = { "tab", "qnum", "snum", "qcovhsp", "scovhsp", "corrected_bitscore" };
	}
	else {
		config.output_format = { "tab", "qseqid", "sseqid", "qcovhsp", "scovhsp", "corrected_bitscore" };
		config.qnum_offset = config.snum_offset = 0;
	}
	if (bi == bj) {
		config.lin_stage1_query = true;
		config.self = true;
		config.query_file.clear();
		config.lin_stage1_combo = false;
	} else {
		config.query_file = { volumes[bi].path };
		config.lin_stage1_query = false;
		config.self = false;
		config.lin_stage1_combo = true;
	}
	config.algo = Config::Algo::DOUBLE_INDEXED;
	config.max_target_seqs_ = INT64_MAX;	
	config.mapany = false;
	config.lin_stage1_target = false;
	config.lowmem_ = 1;
	config.chunk_size = 1024;
	config.database = volumes[bj].path;
	config.db_size = 1000000000;
	config.comp_based_stats = 0;	
	config.output_file = base_dir + std::to_string(bi) + "_" + std::to_string(bj) + ".tsv";
	log_rss();
	Search::run(nullptr, nullptr, nullptr);
}

static string run_block_combos(Job& job, const VolumedFile& volumes) {
	const int64_t BUF_SIZE = 4096;		
	const string base_dir = job.base_dir() + PATH_SEPARATOR + "alignments" + PATH_SEPARATOR, qpath = base_dir + "queue";
	if (job.round() == 0)
		mkdir(job.root_dir() + "accessions");
	mkdir(base_dir);
	Atomic q(qpath);
	const int64_t n = (int64_t)volumes.size();
	const bool accs = job.round() == 0;
	int64_t r;
	atomic<uint64_t> combos_processed(0);
	while (r = q.fetch_add(), r < combos(n)) {
		int64_t bi, bj;
		std::tie(bi, bj) = rank_to_combo(r, n);
		job.log("Searching blocks. Rank=%lli/%lli Blocks=%lli,%lli", r + 1, combos(n), bi, bj);
		if (accs && bi == bj) {
			const string name = job.root_dir() + "accessions" + PATH_SEPARATOR + std::to_string(bi) + ".txt";
			ofstream acc_out(name);
			if (!acc_out.good())
				throw runtime_error("Error opening file " + name);
			unique_ptr<SequenceFile> in(SequenceFile::auto_create({ volumes[bi].path }));
			string id;
			vector<Letter> seq;
			while (in->read_seq(seq, id, nullptr))
				acc_out << Util::Seq::seqid(id.c_str()) << endl;
			in->close();
		}
		run_block_combo(job, volumes, bi, bj, base_dir);
		combos_processed.fetch_add(1, std::memory_order_relaxed);
	}
	Atomic finished(base_dir + "finished");
	finished.fetch_add(combos_processed);
	finished.await(combos(n));
	Atomic concat_lock(base_dir + "concat_lock");
	Atomic concat_done(base_dir + "concat_done");
	if (concat_lock.fetch_add() == 0) {
		const string aln_path = job.base_dir() + PATH_SEPARATOR + "alignments.tsv";
		job.log("Concatenating alignment files to %s", aln_path.c_str());
		{
			ofstream out(aln_path, std::ios::binary | std::ios::trunc);
			if (!out.good())
				throw runtime_error("Error opening file " + aln_path);
			for (int64_t i = 0; i < n; ++i) {
				for (int64_t j = i; j < n; ++j) {
					const string src = base_dir + std::to_string(i) + "_" + std::to_string(j) + ".tsv";
					std::ifstream in(src, std::ios::binary);
					if (!in.good())
						throw runtime_error("Error opening file " + src);
					if (in.peek() == std::ifstream::traits_type::eof())
						continue;
					out << in.rdbuf();
					if (!out) throw runtime_error("Error writing " + aln_path);
				}
			}
		}
		/*const string acc_path = job.root_dir() + PATH_SEPARATOR + "accessions.txt";
		if (job.round() == 0) {
			job.log("Concatenating accession files to %s", acc_path.c_str());
			ofstream out(acc_path);
			if (!out.good())
				throw runtime_error("Error opening file " + acc_path);
			for (int64_t i = 0; i < n; ++i) {
				const string src = job.root_dir() + "accessions" + PATH_SEPARATOR + std::to_string(i) + ".txt";
				std::ifstream in(src);
				if (in.good())
					out << in.rdbuf();
				if (!out) throw runtime_error("Error writing " + acc_path);
			}
		}*/
		const string acc_path = job.root_dir() + PATH_SEPARATOR + "oids.txt";
		if (job.round() == 0) {
			ofstream oid_out(acc_path);
			job.log("Writing oid file");
			for (OId i = 0; i <= volumes.max_oid(); ++i)
				oid_out << i << endl;
		}
		job.log("Running greedy vertex cover");
		config.database = job.round() == 0 ? acc_path : job.root_dir() + "round" + std::to_string(job.round() - 1) + PATH_SEPARATOR + "rep_ids";
		config.edges = aln_path;
		config.edge_format.clear();
		config.output_file = job.base_dir() + PATH_SEPARATOR + "clusters.tsv";
		greedy_vertex_cover();
		if (!job.last_round()) {
			job.log("Writing representative accessions");
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
		concat_done.fetch_add();
	}
	else
		concat_done.await(1);
	return get_reps(job, volumes);
}

string round(Job& job, const VolumedFile& volumes) {
	if (config.mutual_cover.present()) {
		config.min_length_ratio = config.sensitivity < Sensitivity::LINCLUST_40 ?
			std::min(config.mutual_cover.get_present() / 100 + 0.05, 1.0)
			: config.mutual_cover.get_present() / 100 - 0.05;
	}
	job.log("Starting round %i sensitivity %s", job.round(), to_string(config.sensitivity).c_str());
	job.set_round(volumes.sparse_records());
	return run_block_combos(job, volumes);
}

void multinode() {
	if (config.output_file.empty())
		throw runtime_error("Option missing: output file (--out/-o)");
	const string output_file = config.output_file;
	config.file_buffer_size = 64 * 1024;
	
	VolumedFile volumes(config.database.get_present());
	Job job(volumes.max_oid(), volumes.size());
	if (job.worker_id() == 0) {
		if (config.mutual_cover.present())
			job.log("Bi-directional coverage = %f", config.mutual_cover.get_present());
		else
			job.log("Uni-directional coverage = %f", config.member_cover.get(80));
		job.log("Approx. id = %f", config.approx_min_id.get(0));
		job.log("#Volumes = %lli", volumes.size());
		job.log("#Sequences = %lli", volumes.sparse_records());
	}
	if (config.mutual_cover.present()) {
		config.query_or_target_cover = 0;
		config.query_cover = config.mutual_cover.get_present();
		config.subject_cover = config.mutual_cover.get_present();
	}
	else {
		config.query_or_target_cover = config.member_cover.get(80);
		config.query_cover = 0;
		config.subject_cover = 0;
	}
#ifdef WIN32
	_setmaxstdio(8192);
#endif
	vector<string> steps = Cluster::cluster_steps(config.approx_min_id.get(0), true);
	string reps;
	job.set_round_count((int)steps.size());
	const double evalue_cutoff = config.max_evalue,
		target_approx_id = config.approx_min_id.get_present();

	Atomic lock(job.root_dir() + "startup_lock"), done(job.root_dir() + "startup_done");
	if (lock.fetch_add() == 0) {
		const string index_dir = job.root_dir() + "index" + PATH_SEPARATOR;
		mkdir(index_dir);
		for (std::vector<Volume>::const_iterator i = volumes.begin(); i != volumes.end(); ++i) {
			job.log("Indexing volume %td/%zu", i - volumes.begin(), volumes.size());
			FastaFile::index(i->path, index_dir + std::to_string(i - volumes.begin()));
		}
		done.fetch_add();
	}
	else
		done.await(1);

	for (size_t i = 0; i < steps.size(); ++i) {
		config.sensitivity = from_string<Sensitivity>(rstrip(steps[i], "_lin"));
		const vector<string> round_approx_id = config.round_approx_id.empty() ? Cluster::default_round_approx_id(job.round_count()) : config.round_approx_id;
		config.approx_min_id = std::max(target_approx_id, Cluster::round_value(round_approx_id, "--round-approx-id", job.round(), job.round_count()));
		config.max_evalue = i == steps.size() - 1 ? evalue_cutoff : std::min(evalue_cutoff, CASCADED_ROUND_MAX_EVALUE);
		reps = round(job, i == 0 ? volumes : VolumedFile(reps));
		if (i < steps.size() - 1)
			job.next_round();
	}
	Atomic output_lock(job.root_dir() + PATH_SEPARATOR + "output_lock");
	config.output_file = output_file;
	if (output_lock.fetch_add() == 0)
		merge(job, volumes);
}