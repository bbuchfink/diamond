#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>

using std::string;

struct Config
{
	string	input_ref_file;
	unsigned	threads_;
	string	database;
	string	query_file;
	unsigned	merge_seq_treshold;
	unsigned	hit_cap;
	int		min_ungapped_raw_score;
	unsigned shapes;
	unsigned	index_mode;
	unsigned long	max_alignments;
	string	match_file1;
	string	match_file2;
	int		padding;
	unsigned	output_threads;
	unsigned compression;
	unsigned		lowmem;
	double	chunk_size;
	unsigned min_identities;
	unsigned min_identities2;
	int		xdrop;
	unsigned window;
	int		min_hit_score;
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
	string	output_format;
	string	output_file;
	bool		forwardonly;
	unsigned fetch_size;
	bool		single_domain;
	unsigned long	db_size;
	double	query_cover;

	bool		mode_sensitive;
	unsigned	verbosity;
	bool no_auto_append;
	unsigned local_align_mode;

	typedef enum { makedb = 0, blastp = 1, blastx = 2, view = 3, help = 4, version = 5 } Command;
	unsigned	command;

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

	void set_chunk_size(double x);

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
			if (len <= 35)
				return 5;
			else if (len <= 55)
				return 16;
			else
				return 32;
		}
		else
			return padding;
	}

	bool mem_buffered() const
	{
		return tmpdir == "/dev/shm";
	}

	template<typename _t>
	static void set_option(_t& option, _t value)
	{
		if (option == 0)
			option = value;
	}

};

extern Config config;

#endif