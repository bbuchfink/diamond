/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include <string>
#include <vector>
#include <stdint.h>

enum class Sensitivity { FAST = 0, MID_SENSITIVE = 1, SENSITIVE = 2, MORE_SENSITIVE = 3, VERY_SENSITIVE = 4, ULTRA_SENSITIVE = 5 };
enum class TracebackMode { NONE = 0, SCORE_ONLY = 1, STAT = 2, VECTOR = 3, SCORE_BUFFER = 4 };

struct Config
{

	using string = std::string;
	using string_vector = std::vector<std::string>;

	string_vector input_ref_file;
	unsigned	threads_;
	string	database;
	string	query_file;
	unsigned	merge_seq_treshold;
	unsigned	hit_cap;
	unsigned shapes;
	unsigned	index_mode;
	uint64_t	max_alignments;
	string	match_file1;
	string	match_file2;
	int		padding;
	unsigned	output_threads;
	unsigned compression;
	unsigned		lowmem;
	double	chunk_size;
	unsigned min_identities;
	unsigned min_identities2;
	double ungapped_xdrop;
	int		raw_ungapped_xdrop;
	unsigned	min_compressed_identities;
	int		min_seed_score;
	unsigned	seed_signatures;
	double	min_bit_score;
	unsigned	run_len;
	double	max_seed_freq;
	string	tmpdir;
	string	parallel_tmpdir;
	bool		long_mode;
	int		gapped_xdrop;
	double	max_evalue;
	string	kegg_file;
	int		gap_open;
	int		gap_extend;
	string	matrix;
	bool debug_log, verbose, quiet;
	bool		salltitles;
	int		reward;
	int		penalty;
	string	db_type;
	double	min_id;
	unsigned	compress_temp;
	double	toppercent;
	string	daa_file;
	string_vector	output_format;
	string	output_file;
	bool		forwardonly;
	unsigned fetch_size;
	uint64_t	db_size;
	double	query_cover;

	bool		mode_sensitive;
	unsigned	verbosity;
	bool no_auto_append;
	unsigned local_align_mode;
	bool extend_all;
	bool slow_search;
	string_vector seq_no;
	double rank_ratio;
	double rank_factor;
	bool ht_mode;
	bool old_freq;
	double freq_sd;
	unsigned target_fetch_size;
	bool mode_more_sensitive;
	string matrix_file;
	double lambda, K;
	string_vector shape_mask;
	unsigned query_gencode;
	string unaligned;
	double space_penalty;
	bool new_prefilter;
	bool reverse;
	unsigned comp_based_stats;
	int neighborhood_score;
	unsigned seed_weight;
	int report_unaligned;
	double subject_cover;
	bool mode_very_sensitive;
	unsigned max_hsps;
	bool no_self_hits;
	unsigned id_left, id_right, id_n;
	int bmatch, bmismatch, bcutoff;
	unsigned query_bins;
	uint64_t n_ants;
	double rho;
	double p_best;
	double d_exp, d_new;
	double score_estimate_factor;
	int diag_min_estimate;
	string qfilt, sfilt;
	double path_cutoff;
	bool use_smith_waterman;
	string prot_accession2taxid;
	int superblock;
	unsigned max_cells;
	int masking;
	bool benchmark_ranking;
	bool log_query;
	bool log_subject;
	unsigned threads_align;
	double score_ratio;
	bool small_query;
	bool hashed_seeds;
	string nodesdmp;
	bool sallseqid;
	string query_strands;
	bool xml_blord_format;
	int frame_shift;
	bool query_range_culling;
	double query_range_cover;
	double transcript_len_estimate;
	string family_counts_file;
	string taxonlist;
	bool radix_cluster_buffered;
	unsigned join_split_size;
	unsigned join_split_key_len;
	unsigned radix_bits;
	double join_ht_factor;
	bool sort_join;
	bool simple_freq;
	double freq_treshold;
	bool use_lazy_dict;
	string aligned_file;
	int filter_locus;
	bool use_dataset_field;
	bool store_query_quality;
	string invocation;
	unsigned swipe_chunk_size;
	unsigned query_parallel_limit;
	bool long_reads;
	bool output_header;
	string alfmt;
	string unfmt;
	string namesdmp;
	bool hardmasked;
	int cbs_window;
	double tantan_r;
	double tantan_minMaskProb;
	bool no_unlink;
	bool no_dict;
	int stop_match_score;
	int tantan_maxRepeatOffset;
	bool tantan_ungapped;
	string taxon_exclude;
	bool swipe_all;
	uint64_t taxon_k;
	uint64_t upgma_edge_limit;
	string tree_file;
	string upgma_dist;
	string upgma_input;
	bool log_extend;
	int chaining_maxgap;
	string family_map;
	string family_map_query;
	size_t chaining_range_cover;
	bool no_swipe_realign;
	bool cut_bar;
	bool bootstrap;
	size_t chaining_maxnodes;
	int cutoff_score_8bit;
	double inner_culling_overlap;
	double min_band_overlap;
	int min_realign_overhang;
	int ungapped_window;
	int gapped_filter_diag_score;
	double gapped_filter_evalue;
	int gapped_filter_window;
	bool output_hits;
	double ungapped_evalue;
	bool no_logfile;
	bool no_heartbeat;
	int band_bin;
	int col_bin;
	size_t file_buffer_size;
	bool self;
	size_t trace_pt_fetch_size;
	uint32_t tile_size;
	double short_query_ungapped_bitscore;
	int short_query_max_len;
	double gapped_filter_evalue1;
	size_t ext_chunk_size;
	double ext_min_yield;
	string ext;
	int full_sw_len;
	double relaxed_evalue_factor;
	string type;
	bool raw;
	bool mode_ultra_sensitive;
	double chaining_len_cap;
	size_t chaining_min_nodes;
	bool fast_tsv;
	unsigned target_parallel_verbosity;
	double memory_limit;
	size_t global_ranking_targets;
	bool mode_mid_sensitive;

	Sensitivity sensitivity;
	TracebackMode traceback_mode;

	bool multiprocessing;
	bool mp_init;

	enum {
		makedb = 0, blastp = 1, blastx = 2, view = 3, help = 4, version = 5, getseq = 6, benchmark = 7, random_seqs = 8, compare = 9, sort = 10, roc = 11, db_stat = 12, model_sim = 13,
		match_file_stat = 14, model_seqs = 15, opt = 16, mask = 17, fastq2fasta = 18, dbinfo = 19, test_extra = 20, test_io = 21, db_annot_stats = 22, read_sim = 23, info = 24, seed_stat = 25,
		smith_waterman = 26, cluster = 27, translate = 28, filter_blasttab = 29, show_cbs = 30, simulate_seqs = 31, split = 32, upgma = 33, upgma_mc = 34, regression_test = 35,
		reverse_seqs = 36, compute_medoids = 37, mutate = 38, merge_tsv = 39
	};
	unsigned	command;

	enum { double_indexed = 0, query_indexed = 1, subject_indexed = 2 };
	int algo;

	string cluster_algo;
	double cluster_mcl_inflation;
	uint32_t cluster_mcl_expansion;
	double cluster_mcl_sparsity_switch;
	string	cluster_similarity;

	enum { query_parallel = 0, target_parallel = 1 };
	unsigned load_balancing;

	Config() {}
	Config(int argc, const char **argv, bool check_io = true);

	inline unsigned get_run_len(unsigned length)
	{
		if (run_len == 0) {
			if (length < 30)
				return 1;
			else if (length < 100)
				return 20;
			else
				return 40;
		}
		else
			return run_len;
	}

	inline bool output_range(unsigned n_target_seq, int score, int top_score)
	{
		if (toppercent < 100)
			return (1.0 - (double)score / top_score) * 100 <= toppercent;
		else
			return n_target_seq < max_alignments;
	}

	void set_sens(Sensitivity sens);

	bool mem_buffered() const { return tmpdir == "/dev/shm"; }

  	template<typename _t>
	static void set_option(_t& option, _t value, _t def = 0) { if (option == def) option = value; }
};

void print_warnings();
extern Config config;