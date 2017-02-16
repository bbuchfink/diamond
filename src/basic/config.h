/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
	string	seg;
	bool debug_log, verbose, quiet;
	bool		have_ssse3;
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
	double rank_factor;
	double rank_ratio;
	bool ht_mode;
	bool old_freq;
	double freq_sd;
	bool query_parallel;
	unsigned target_fetch_size;
	bool mode_more_sensitive;
	string matrix_file;
	double lambda, K;
	vector<string> shape_mask;
	unsigned seed_anchor;
	unsigned query_gencode;
	string unaligned;
	double space_penalty, raw_space_penalty;
	double min_diag_score;
	int min_diag_raw_score;
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
	bool greedy;
	string qfilt;

	enum {
		makedb = 0, blastp = 1, blastx = 2, view = 3, help = 4, version = 5, getseq = 6, benchmark = 7, random_seqs = 8, compare = 9, sort = 10, roc = 11, db_stat = 12, model_sim = 13,
		match_file_stat = 14, model_seqs = 15, opt = 16
	};
	unsigned	command;

	enum { double_indexed = 0, subject_indexed = 1 };
	unsigned algo;

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