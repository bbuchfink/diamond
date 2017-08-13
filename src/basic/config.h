/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>
#include <vector>
#include <stdint.h>

using std::string;
using std::vector;

struct Config
{
	string	input_ref_file;
	unsigned	threads_;
	string	database;
	string	query_file;
	unsigned	merge_seq_treshold;
	unsigned	hit_cap;
	double min_ungapped_score;
	int		min_ungapped_raw_score;
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
	unsigned window;
	double		min_hit_score;
	int min_hit_raw_score;
	int		hit_band;
	unsigned	min_compressed_identities;
	int		min_seed_score;
	unsigned	seed_signatures;
	double	min_bit_score;
	unsigned	run_len;
	bool		alignment_traceback;
	double	max_seed_freq;
	string	tmpdir;
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
	vector<string>	output_format;
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
	vector<string> seq_no;
	double rank_ratio;
	double rank_ratio2;
	bool ht_mode;
	bool old_freq;
	double freq_sd;
	unsigned target_fetch_size;
	bool mode_more_sensitive;
	string matrix_file;
	double lambda, K;
	vector<string> shape_mask;
	unsigned seed_anchor;
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
	double query_overlap_culling;
	string query_strands;

	enum {
		makedb = 0, blastp = 1, blastx = 2, view = 3, help = 4, version = 5, getseq = 6, benchmark = 7, random_seqs = 8, compare = 9, sort = 10, roc = 11, db_stat = 12, model_sim = 13,
		match_file_stat = 14, model_seqs = 15, opt = 16, mask = 17, fastq2fasta=18, dbinfo=19
	};
	unsigned	command;

	enum { double_indexed = 0, query_indexed = 1, subject_indexed = 2 };
	int algo;

	enum { query_parallel = 0, target_parallel = 1 };
	unsigned load_balancing;

	enum {
		swipe = 0, greedy = 1, floating_xdrop = 4, more_greedy = 2, most_greedy=3,
	};
	int ext;

	Config() {}
	Config(int argc, const char **argv);

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

	/*unsigned read_padding(size_t len)
	{
		if (padding == 0) {
			if (len <= 255)
				return 10;
			else
				return 32;
		}
		else
			return padding;
	}*/

	unsigned read_padding(size_t len)
	{
		if (padding == 0) {
			if (mode_very_sensitive)
				return 60;
			else if (len <= 35)
				return 5;
			else if (len <= 55)
				return 16;
			else
				return 32;
		}
		else
			return padding;
	}

	bool mem_buffered() const { return tmpdir == "/dev/shm"; }

  	template<typename _t>
	static void set_option(_t& option, _t value) { if (option == 0) option = value; }
};

extern Config config;

#endif