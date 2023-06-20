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
#include <map>
#include "../util/enum.h"
#include "../util/options/option.h"
#include "value.h"

enum class Sensitivity { FASTER = -1, FAST = 0, DEFAULT = 1, MID_SENSITIVE = 2, SENSITIVE = 3, MORE_SENSITIVE = 4, VERY_SENSITIVE = 5, ULTRA_SENSITIVE = 6};
enum class Compressor;

template<> struct EnumTraits<Sensitivity> {
	static const EMap<Sensitivity> to_string;
	static const SEMap<Sensitivity> from_string;
};

template<> struct EnumTraits<SequenceType> {
    static const EMap<SequenceType> to_string;
    static const SEMap<SequenceType> from_string;
};

enum class GraphAlgo { GREEDY_VERTEX_COVER, LEN_SORTED };

template<> struct EnumTraits<GraphAlgo> {
	static const SEMap<GraphAlgo> from_string;
};

struct CommandLineParser;

constexpr int64_t DEFAULT_MAX_TARGET_SEQS = 25;

struct Config
{

	using string = std::string;
	using string_vector = std::vector<std::string>;

	string_vector input_ref_file;
	int	threads_;
	Option<string> database;
	string_vector query_file;
	unsigned	merge_seq_treshold;
	unsigned shapes;
	Option<int64_t> max_target_seqs_;
	string	match_file1;
	string	match_file2;
	int		padding;
	unsigned	output_threads;
	string compression;
	unsigned		lowmem_;
	double	chunk_size;
	unsigned min_identities_;
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
	double gapped_xdrop;
	double	max_evalue;
	string	kegg_file;
	int		gap_open;
	int		gap_extend;
    int mismatch_penalty;
    int match_reward;
	string	matrix;
	bool debug_log, verbose, quiet;
	bool		salltitles;
	int		reward;
	int		penalty;
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
	double query_or_target_cover;
	bool		mode_sensitive;
	unsigned	verbosity;
	bool no_auto_append;
	string_vector seq_no;
	double rank_factor;
	double rank_ratio;
	double freq_sd_;
	unsigned target_fetch_size;
	bool mode_more_sensitive;
	string matrix_file;
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
	unsigned query_bins_;
	uint64_t n_ants;
	double rho;
	double p_best;
	double d_exp, d_new;
	double score_estimate_factor;
	int diag_min_estimate;
	double path_cutoff;
	bool use_smith_waterman;
	string prot_accession2taxid;
	int superblock;
	unsigned max_cells;
	Option<string> masking_;
	bool log_query;
	bool log_subject;
	unsigned threads_align;
	double score_ratio;
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
	string aligned_file;
	bool use_dataset_field;
	bool store_query_quality;
	string invocation;
	unsigned swipe_chunk_size;
	unsigned query_parallel_limit;
	bool long_reads;
	Option<string_vector> output_header;
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
	double gapped_filter_diag_bit_score;
	double gapped_filter_evalue_;
	int gapped_filter_window;
	bool output_hits;
	double ungapped_evalue_;
	bool no_logfile;
	int band_bin;
	int col_bin;
	size_t file_buffer_size;
	bool self;
	int64_t trace_pt_fetch_size;
	uint32_t tile_size;
	double short_query_ungapped_bitscore;
	int short_query_max_len;
	double gapped_filter_evalue1;
	size_t ext_chunk_size;
	double ext_min_yield;
	string ext_;
	int full_sw_len;
	double relaxed_evalue_factor;
	string type;
	bool raw;
	bool mode_ultra_sensitive;
	double chaining_len_cap;
	size_t chaining_min_nodes;
	bool fast_tsv;
	unsigned target_parallel_verbosity;
	int64_t global_ranking_targets;
	bool mode_mid_sensitive;
	bool no_ranking;
	bool query_memory;
	size_t memory_intervals;
	size_t seedhit_density;
	size_t chunk_size_multiplier;
	double ranking_score_drop_factor;
	double ranking_cutoff_bitscore;
	int left_most_interval;
	bool no_forward_fp;
	bool no_ref_masking;
	string roc_file;
	bool target_bias;
	bool check_multi_target;
	bool output_fp;
	int family_cap;
	int cbs_matrix_scale;
	size_t query_count;
	double cbs_err_tolerance;
	int cbs_it_limit;
	double query_match_distance_threshold;
	double length_ratio_threshold;
	bool hash_join_swap;
	bool target_indexed;
	size_t deque_bucket_size;
	bool mode_fast;
	double log_evalue_scale;
	double ungapped_evalue_short_;
	int64_t max_swipe_dp;
	std::string seqidlist;
	bool skip_missing_seqids;
	Option<string_vector> iterate;
	bool ignore_warnings;
	bool short_seqids;
	bool no_reextend;
	double seed_cut_;
	bool no_reorder;
	string file1;
	string file2;
	size_t key2;
	string motif_mask_file;
	string motif_masking;
	bool freq_masking;
	Loc max_motif_len;
	double chaining_stacked_hsp_ratio;
	Option<double> cluster_threshold;
	Option<string> memory_limit;
	int64_t swipe_task_size;
	Loc minimizer_window_;
	bool lin_stage1;
	int64_t min_task_trace_pts;
	Loc sketch_size;
	string soft_masking;
	string oid_list;
	int64_t bootstrap_block;
	int64_t centroid_factor;
	int timeout;
	string resume;
	int64_t target_hard_cap;
	bool mapany;
	Option<string> clustering;
	Option<string> neighbors;
	double reassign_overlap;
	double reassign_ratio;
	int64_t reassign_max;
	bool add_self_aln;
	string centroid_out;
	string unaligned_targets;
	Option<double> approx_min_id;
	bool mode_faster;
	double member_cover;
	bool weighted_gvc;
	bool kmer_ranking;
	bool hamming_ext;
	double diag_filter_id;
	double diag_filter_cov;
	bool strict_gvc;
	bool mmseqs_compat;
	string edge_format;
	bool no_block_size_limit;
	string edges;
	bool mp_self;
	bool approx_backtrace;
	bool prefix_scan;
	double narrow_band_cov;
	double narrow_band_factor;
	Loc anchor_window;
	double anchor_score;
	bool classic_band;
	bool no_8bit_extension;
	bool anchored_swipe;
	bool no_chaining_merge_hsps;
	bool recluster_bd;
	bool pipeline_short;
	string graph_algo;
	bool linsearch;
	int64_t tsv_read_size;
    int zdrop;
	bool heartbeat;
	bool no_parse_seqids;
	bool sam_qlen_field;

    SequenceType dbtype;

	Sensitivity sensitivity;

	bool multiprocessing;
	bool mp_init;
	bool mp_recover;
	int mp_query_chunk;

	enum {
		makedb = 0, blastp = 1, blastx = 2, view = 3, help = 4, version = 5, getseq = 6, benchmark = 7, random_seqs = 8, compare = 9, sort = 10, roc = 11, db_stat = 12, model_sim = 13,
		match_file_stat = 14, model_seqs = 15, opt = 16, mask = 17, fastq2fasta = 18, dbinfo = 19, test_extra = 20, test_io = 21, db_annot_stats = 22, read_sim = 23, info = 24, seed_stat = 25,
		smith_waterman = 26, cluster = 27, translate = 28, filter_blasttab = 29, show_cbs = 30, simulate_seqs = 31, split = 32, upgma = 33, upgma_mc = 34, regression_test = 35,
		reverse_seqs = 36, compute_medoids = 37, mutate = 38, rocid = 40, makeidx = 41, find_shapes, prep_db, composition, JOIN, HASH_SEQS, LIST_SEEDS, CLUSTER_REALIGN,
		GREEDY_VERTEX_COVER, INDEX_FASTA, FETCH_SEQ, CLUSTER_REASSIGN, blastn, RECLUSTER, LENGTH_SORT, MERGE_DAA, DEEPCLUST, LINCLUST, WORD_COUNT, CUT, MODEL_SEQS
	};


	unsigned	command;

	enum class Algo { AUTO = -1, DOUBLE_INDEXED = 0, QUERY_INDEXED = 1, CTG_SEED };
	Algo algo;

	Option<string> cluster_algo;
	Option<string> cluster_similarity;
	string cluster_graph_file;
	bool cluster_restart;

	size_t max_size_set;
	bool external;
	string_vector cluster_steps;
	double cluster_mcl_inflation;
	uint32_t cluster_mcl_expansion;
	double cluster_mcl_sparsity_switch;
	uint32_t cluster_mcl_chunk_size;
	uint32_t cluster_mcl_max_iter;
	bool cluster_mcl_stats;
	bool cluster_mcl_nonsymmetric;

	enum { query_parallel = 0, target_parallel = 1 };
	unsigned load_balancing;

	Config() {}
	Config(int argc, const char **argv, bool check_io, CommandLineParser& parser);

	Loc min_orf_len(Loc length) const {
		if (run_len == 0) {
			if (length < 30 || frame_shift != 0)
				return 1;
			else if (length < 100)
				return 20;
			else
				return 40;
		}
		else
			return run_len;
	}

	inline bool output_range(unsigned n_target_seq, int score, int top_score, const int64_t max_target_seqs)
	{
		if (toppercent < 100)
			return (1.0 - (double)score / top_score) * 100 <= toppercent;
		else
			return n_target_seq < max_target_seqs;
	}

	int64_t block_size() const {
		return (int64_t)(chunk_size * 1e9);
	}

	void set_sens(Sensitivity sens);
	std::string single_query_file() const;

	bool mem_buffered() const { return tmpdir == "/dev/shm"; }
	Compressor compressor() const;

  	template<typename _t>
	static void set_option(_t& option, _t value, _t def = 0) { if (option == def) option = value; }
	template<typename _t>
	static void set_option(_t& option, _t value, _t def, _t alt) { if (value != def) option = value; else option = alt; }
};

void print_warnings();
extern Config config;

template<typename _t>
_t top_cutoff_score(_t top_score) {
	return _t((1.0 - config.toppercent / 100.0)*top_score);
}

template<> struct EnumTraits<Config::Algo> {
	static const EMap<Config::Algo> to_string;
	static const SEMap<Config::Algo> from_string;
};

extern const char* const DEFAULT_MEMORY_LIMIT;

std::pair<double, int> block_size(int64_t memory_limit, Sensitivity s, bool lin);
