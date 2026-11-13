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
#include <ctime>
#include <cstring>
#include <exception>
#include <iostream>
#include <algorithm>
#include <string>
#include "basic/config.h"
#include "volume.h"
#include "multinode.h"
#include "data/sequence_file.h"
#include "data/block/block.h"
#include "cluster/cluster.h"
#include "tools/tools.h"
#include "util/log_stream.h"
#include "util/io/file.h"
#include "util/command_line_parser.h"
#include "util/parallel/multiprocessing.h"

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

static const std::chrono::minutes HEARTBEAT_INTERVAL(15);
static const uint32_t HEARTBEAT_TIMEOUT_S = 20 * 60;

Heartbeat::Heartbeat(const string& file_name, Job& job) :
	file_(file_name),
	job_(job),
	worker_id_(job.worker_id()),
	stop_(false),
	failed_(false),
	uncaught_exceptions_(std::uncaught_exceptions()),
	thread_(&Heartbeat::run, this)
{}

Heartbeat::~Heartbeat() {
	{
		std::lock_guard<std::mutex> lock(mtx_);
		stop_ = true;
	}
	cv_.notify_one();
	thread_.join();
	if (std::uncaught_exceptions() > uncaught_exceptions_ || failed_.load(std::memory_order_acquire))
		return;
	try {
		const uint32_t t = DONE;
		file_.write_at(worker_id_ * (int64_t)sizeof(t), (const char*)&t, sizeof(t));
	}
	catch (const std::exception& e) {
		std::cerr << "Warning: failed to write heartbeat done marker: " << e.what() << endl;
	}
}

void Heartbeat::run() {
	std::unique_lock<std::mutex> lock(mtx_);
	do {
		const uint32_t t = (uint32_t)std::time(nullptr);
		try {
			if (file_.write_at(worker_id_ * (int64_t)sizeof(t), (const char*)&t, sizeof(t)) != (int64_t)sizeof(t))
				throw runtime_error("Error writing heartbeat file " + file_.file_name());
			update_stale_workers(t);
		}
		catch (...) {
			exception_ = std::current_exception();
			failed_.store(true, std::memory_order_release);
			return;
		}
	} while (!cv_.wait_for(lock, HEARTBEAT_INTERVAL, [this] { return stop_; }));
}

void Heartbeat::update_stale_workers(uint32_t now) {
	const string data = file_.read_all();
	vector<int64_t> stale;
	for (size_t i = 0; i < data.size() / sizeof(uint32_t); ++i) {
		if ((int64_t)i == worker_id_)
			continue;
		uint32_t t;
		memcpy(&t, data.data() + i * sizeof(uint32_t), sizeof(t));
		// A zero slot belongs to a worker that has not written its first heartbeat yet.
		if (t != 0 && t != DONE && t < now && now - t > HEARTBEAT_TIMEOUT_S)
			stale.push_back((int64_t)i);
	}
	job_.set_stale_workers(std::move(stale));
}

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
	check_heartbeat();
	size_t n = message.size();
	while (n > 0 && (message[n - 1] == '\n' || message[n - 1] == '\r'))
		--n;
	const string s = log_prefix() + message.substr(0, n) + LOG_EOL;
	*message_stream << s;
	log_file_->push(s);
}

bool Job::stopped() {
	check_heartbeat();
	if (stopped_)
		return true;
	if (!heartbeat_)
		return false;
	const vector<int64_t> stale = stale_workers();
	if (!stale.empty()) {
		log("Stale workers detected (no heartbeat for more than 20 minutes): %s, terminating...", join(",", stale.begin(), stale.end()).c_str());
		File(stop_file(), "wb").close();
		stopped_ = true;
	}
	else if (file_exists(stop_file())) {
		log("Stop file %s found, terminating...", stop_file().c_str());
		stopped_ = true;
	}
	return stopped_;
}

void Job::reset_workers() {
	for (const string& file : { stop_file(), base_dir_ + PATH_SEPARATOR + "heartbeat" }) {
		if (!file_exists(file))
			continue;
		if (std::remove(file.c_str()) != 0)
			throw runtime_error("Error removing file " + file);
		log("Recover: removed %s", file.c_str());
	}
	const string worker_id_file = base_dir_ + PATH_SEPARATOR + "worker_id";
	if (file_exists(worker_id_file)) {
		Atomic worker_id(worker_id_file, *this);
		const int64_t n = worker_id.get();
		worker_id.reset();
		log("Recover: reset worker id from %lli to 0", (long long)n);
	}
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

static int64_t finished_count(Job& job, const string& finished_file) {
	if (!file_exists(finished_file))
		return 0;
	AtomicList finished(finished_file, job);
	const vector<int64_t> v = finished.values();
	return (int64_t)std::unordered_set<int64_t>(v.begin(), v.end()).size();
}

bool queue_finished(Job& job, const string& finished_file, int64_t n) {
	return finished_count(job, finished_file) >= n;
}

static void recover_reset(Job& job, const char* step, const string& file_name) {
	if (!file_exists(file_name))
		return;
	Atomic a(file_name, job);
	const int64_t v = a.get();
	a.reset();
	job.log("Recover: %s: reset %s from %lli to 0", step, file_name.c_str(), (long long)v);
}

bool recover_step(Job& job, const char* step, const string& lock_file, const string& done_file, const string& queue_file) {
	if (file_exists(done_file) && Atomic(done_file, job).get() >= 1) {
		job.log("Recover: %s completed", step);
		return true;
	}
	if ((!lock_file.empty() && file_exists(lock_file)) || (!queue_file.empty() && file_exists(queue_file))) {
		job.log("Recover: %s interrupted", step);
		if (!lock_file.empty())
			recover_reset(job, step, lock_file);
		if (!queue_file.empty())
			recover_reset(job, step, queue_file);
	}
	else
		job.log("Recover: %s not begun", step);
	job.stop();
	return false;
}

bool recover_queue(Job& job, const char* step, const string& queue_file, const string& finished_file, int64_t n) {
	const int64_t finished = finished_count(job, finished_file);
	if (finished >= n) {
		job.log("Recover: %s completed", step);
		return true;
	}
	if (file_exists(queue_file)) {
		job.log("Recover: %s interrupted, %lli/%lli items finished", step, (long long)finished, (long long)n);
		recover_reset(job, step, queue_file);
	}
	else
		job.log("Recover: %s not begun", step);
	job.stop();
	return false;
}

static string seed_hit_table_file(const string& base_dir, int64_t member_superblock, int64_t rep_superblock) {
	return base_dir + "seed_hits_" + std::to_string(member_superblock) + "_" + std::to_string(rep_superblock);
}

/* Writes the seed hit table to a temporary file which is then renamed, so that a file under the final
   name is always complete. */
static void write_seed_hit_table(const vector<BitVector>& table, const string& file_name) {
	const string tmp = file_name + ".tmp";
	{
		ofstream out;
		out.exceptions(std::ios::failbit | std::ios::badbit);
		out.open(tmp, std::ios::binary | std::ios::trunc);
		const uint64_t n = table.size();
		out.write((const char*)&n, sizeof(n));
		for (const BitVector& v : table)
			v.write(out);
		out.close();
	}
	if (std::rename(tmp.c_str(), file_name.c_str()) != 0)
		throw runtime_error("Error renaming file " + tmp + " to " + file_name);
}

static void read_seed_hit_table(vector<BitVector>& table, const string& file_name) {
	ifstream in;
	in.exceptions(std::ios::failbit | std::ios::badbit);
	in.open(file_name, std::ios::binary);
	uint64_t n;
	in.read((char*)&n, sizeof(n));
	table.clear();
	table.resize((size_t)n);
	for (BitVector& v : table)
		v.read(in);
}

static void remove_if_exists(const string& file_name) {
	if (file_exists(file_name))
		remove_tmp_file(file_name);
}

static void run_block_combos(Job& job, const VolumedFile& superblocks, const string& base_dir, const string& aln_volumes_path) {
	int64_t member_superblock;
	const bool lin_index = use_lin_index(job);
	const int64_t n = (int64_t)superblocks.size();
	// The seed indices are not needed any more (and may have been removed) once all block combinations are finished.
	if (lin_index && !(config.recover && queue_finished(job, base_dir + "finished", n))) {
		configure_round(job, superblocks.letter_count());
		build_lin_indices(job, superblocks);
		if (!job.goon())
			return;
	}
	if (config.recover && !recover_queue(job, "Block search", base_dir + "queue", base_dir + "finished", n))
		return;
	AtomicList finished(base_dir + "finished", job);
	Atomic q(base_dir + "queue", job, finished);
	while (!config.recover && job.goon() && (member_superblock = q.fetch_add(), member_superblock < n)) {
		unique_ptr<vector<BitVector>> seed_hit_table;
		if (member_superblock > 0)
			seed_hit_table.reset(new vector<BitVector>());
		const shared_ptr<Block> member_block = load_member_block(job, superblocks, member_superblock);
		/* The rep superblocks of a member superblock are processed in order by the single worker that owns
		   the member superblock, since the seed hit table carries over from one rep superblock to the next.
		   The queue only tracks progress, so that an interrupted member superblock resumes at its first
		   unfinished rep superblock. As no other worker touches these files, the owner can reset the queue. */
		const string rep_suffix = "_" + std::to_string(member_superblock);
		AtomicList rep_finished(base_dir + "rep_finished" + rep_suffix, job);
		/* After each rep superblock except the last, the seed hit table is saved as seed_hits_<member>_<rep>
		   before the rep superblock is marked finished. A resumed member superblock restores the table of
		   the last finished rep superblock. A table file one past that means the worker failed after saving
		   it, so that rep superblock's search had completed and it is marked finished here. */
		if (seed_hit_table) {
			const vector<int64_t> done = rep_finished.values();
			int64_t last = done.empty() ? -1 : *std::max_element(done.begin(), done.end());
			if (file_exists(seed_hit_table_file(base_dir, member_superblock, last + 1))) {
				++last;
				rep_finished.push(last);
			}
			if (last >= 0 && last < member_superblock) {
				const string f = seed_hit_table_file(base_dir, member_superblock, last);
				if (!file_exists(f))
					throw runtime_error("Missing seed hit table file " + f);
				read_seed_hit_table(*seed_hit_table, f);
				job.log("Restored seed hit table after rep superblock %lli", (long long)last + 1);
			}
			if (last >= 1)
				remove_if_exists(seed_hit_table_file(base_dir, member_superblock, last - 1));
			remove_if_exists(seed_hit_table_file(base_dir, member_superblock, last + 1) + ".tmp");
		}
		Atomic rep_q(base_dir + "rep_queue" + rep_suffix, job, rep_finished);
		rep_q.reset();
		int64_t rep_superblock;
		while (job.goon() && (rep_superblock = rep_q.fetch_add(), rep_superblock <= member_superblock)) {
			job.log("Searching blocks. Blocks=%lli,%lli", member_superblock + 1, rep_superblock + 1);
			/*if (!seed_hit_table->empty()) {
				for (size_t i = 0; i < seed_hit_table->size(); ++i)
					job.log("Seed hit table paired positions shape %zu: %zu/%zu", i, seed_hit_table->operator[](i).one_count(), seed_hit_table->operator[](i).size());
			}*/
			run_search(job, superblocks, member_superblock, rep_superblock, member_block, base_dir, seed_hit_table);
			if (!job.goon())
				break;
			if (seed_hit_table && rep_superblock < member_superblock)
				write_seed_hit_table(*seed_hit_table, seed_hit_table_file(base_dir, member_superblock, rep_superblock));
			rep_finished.push(rep_superblock);
			if (seed_hit_table && rep_superblock > 0)
				remove_if_exists(seed_hit_table_file(base_dir, member_superblock, rep_superblock - 1));
		}
		if (!job.goon())
			break;
		finished.push(member_superblock);
		rep_q.remove();
		rep_finished.remove();
		if (member_superblock > 0)
			remove_if_exists(seed_hit_table_file(base_dir, member_superblock, member_superblock - 1));
		job.finish_step();
	}
	if (!job.goon())
		return;
	if (!finished.await(n, job))
		return;

	if (config.recover && !recover_step(job, "Alignment volume list", base_dir + "volumes_lock", base_dir + "volumes_done"))
		return;
	Atomic volumes_lock(base_dir + "volumes_lock", job);
	Atomic volumes_done(base_dir + "volumes_done", job);
	if (!config.recover && volumes_lock.fetch_add() == 0) {
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
		volumes_done.await(1, job);
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
		run_search(job, superblocks, -1, -1, nullptr, base_dir, seed_hit_table);
	}	
	if (!job.goon())
		return pair<string, uint64_t>("", 0);
	// Recover mode leaves the removal of files to the next run, which passes through here as well.
	if (!config.recover) {
		superblocks.remove(false, linear && config.mutual_cover.blank(), false);
		if (job.last_round()) {
			if (!config.fasta_index_file.empty())
				remove_tmp_file(config.fasta_index_file);
			superblocks.remove(job.round() > 0, true, false);
		}
	}
	if (config.recover && !recover_step(job, "Greedy vertex cover", base_dir + "gvc_lock", base_dir + "gvc_done"))
		return pair<string, uint64_t>("", 0);
	Atomic gvc_lock(base_dir + "gvc_lock", job);
	Atomic gvc_done(base_dir + "gvc_done", job);
	GVC::Cfg cfg;
	if (!config.recover && gvc_lock.fetch_add() == 0) {
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
	else if (!gvc_done.await(1, job))
		return pair<string, uint64_t>("", 0);
	gvc_lock.close();
	gvc_done.close();
	rmdir(base_dir.c_str());
	if (!job.goon())
		return pair<string, uint64_t>("", 0);
	return get_reps(job, round_minichunks, std::move(cfg.clustering));
}

static string parallel_run_dir(const string& parallel_tmpdir) {
	return parallel_tmpdir + PATH_SEPARATOR + "diamond-tmp-" + Const::version_string + PATH_SEPARATOR;
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
	const bool is_linclust = Cluster::is_linclust(rounds);
	// TODO
	config.hamming_ext = config.approx_min_id.present() ? config.approx_min_id.get_present() >= 50.0 : false;
	//config.freq_masking = true;
	
	if (parallel) {
		config.tmpdir = parallel_run_dir(config.parallel_tmpdir);
		mkdir(config.tmpdir);
	}
	else
		config.tmpdir = create_temp_directory(config.tmpdir, "diamond-tmp-") + PATH_SEPARATOR;	
	Job job;
	if (job.stopped())
		return;
	const string input_vols = config.tmpdir + "input_vols.tsv";
	if (config.recover && !recover_step(job, "Startup", config.tmpdir + "startup_lock", config.tmpdir + "startup_done"))
		return;
	Atomic lock(config.tmpdir + "startup_lock", job), done(config.tmpdir + "startup_done", job);
	if (!config.recover && lock.fetch_add() == 0) {
		ofstream input_vols_file(input_vols);
		input_vols_file << config.database.get_present() << endl;
		input_vols_file.close();
		done.fetch_add();
	}
	else if (!done.await(1, job))
		return;
	config.database = input_vols;
	
	VolumedFile volumes(config.database.get_present());
	if (job.worker_id() == 0) {
		ofstream args_file;
		args_file.exceptions(std::ios::failbit | std::ios::badbit);
		args_file.open(config.tmpdir + "args");
		args_file << config.invocation << endl;
		args_file.close();
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
		const int ccd = Cluster::round_ccd(ccd_arg, i, rounds.size(), linear_round);
		config.connected_component_depth.clear();
		config.connected_component_depth.push_back(std::to_string(ccd));

		const string superblocks = job.round() == 0 ? make_merged_blocks(job, input_minichunks_seqs, job.base_dir() + "input_superblocks" + PATH_SEPARATOR, volumes.letter_count())
			: make_merged_blocks(job, rep_minichunks, job.base_dir() + "input_superblocks" + PATH_SEPARATOR, volumes.letter_count());
		if (!job.goon())
			return;
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
	const string output_lock_file = job.root_dir() + PATH_SEPARATOR + "output_lock", output_done_file = job.root_dir() + PATH_SEPARATOR + "output_done";
	if (config.recover) {
		if (recover_step(job, "Output", output_lock_file, output_done_file))
			job.log("Recover: run is complete");
		return;
	}
	Atomic output_lock(output_lock_file, job), output_done(output_done_file, job);
	config.output_file = output_file;
	if (output_lock.fetch_add() == 0) {
		VolumedFile acc_vols(input_minichunks_accs);
		merge(job, acc_vols, volumes, hdr_format);
		job.log(job.stats());
		output_done.fetch_add();
		output_lock.close();
		output_done.close();
		lock.close();
		done.close();
		if (parallel)
			return;
		job.finish();
		remove_tmp_file(input_vols);
		remove_tmp_file(config.tmpdir + "args");
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

void recover() {
	if (config.parallel_tmpdir.empty())
		throw runtime_error("Option missing: --parallel-tmpdir");
	const string parallel_tmpdir = config.parallel_tmpdir;
	const string args_path = parallel_run_dir(parallel_tmpdir) + "args";
	ifstream args_file(args_path);
	if (!args_file)
		throw runtime_error("Error opening file " + args_path);
	string invocation;
	getline(args_file, invocation);
	args_file.close();
	const vector<string> args = Util::String::tokenize(invocation.c_str(), " \t\r\n");
	vector<const char*> argv;
	for (const string& s : args)
		argv.push_back(s.c_str());
	CommandLineParser parser;
	config = Config((int)argv.size(), argv.data(), true, parser);
	if (config.command != ::Config::cluster && config.command != ::Config::DEEPCLUST && config.command != ::Config::LINCLUST)
		throw runtime_error("Invalid workflow in args file: " + args_path);
	config.parallel_tmpdir = parallel_tmpdir;
	config.recover = true;
	multinode();
}
