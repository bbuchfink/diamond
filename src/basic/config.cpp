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

#include <memory>
#include <cstdlib>
#include <sys/stat.h>
#include <exception>
#include <iomanip>
#include <numeric>
#include "../util/command_line_parser.h"
#include "config.h"
#include "../util/util.h"
#include "../util/log_stream.h"
#include "../basic/value.h"
#include "score_matrix.h"
#include "../util/system.h"
#include "reduction.h"
#include "shape_config.h"
#include "../util/io/temp_file.h"
#include "../basic/match.h"
#include "../cluster/cluster_registry.h"
#include "../basic/translate.h"
#include "../dp/dp.h"
#include "masking.h"
#include "../util/system/system.h"
#include "../util/simd.h"
#include "../util/parallel/multiprocessing.h"

using namespace std;

Config config;

void print_warnings() {
	if (config.sensitivity >= Sensitivity::VERY_SENSITIVE || config.verbosity == 0)
		return;
	const double ram = total_ram();
	unsigned b = 2, c = 4;
	if (ram >= 511) {
		b = 12;
		c = 1;
	} else if (ram >= 255) {
		b = 8;
		c = 1;
	} else if (ram >= 127) {
		b = 5;
		c = 1;
	} else if (ram >= 63) {
		b = 6;
	} else if (ram >= 31) {
		b = 4;
	}
	if ((b > 2 && b > config.chunk_size) || c < config.lowmem) {
		set_color(Color::YELLOW, true);
		cerr << "The host system is detected to have " << (size_t)ram << " GB of RAM. It is recommended to increase the block size for better performance using these parameters: -b" << b;
		if (c != 4)
			cerr << " -c" << c;
		cerr << endl;
		reset_color(true);
	}
}

template<typename _t>
_t set_string_option(const string& s, const string& name, const vector<pair<string, _t>>& values) {
	if (s.empty())
		return (_t)0;
	for (auto& i : values)
		if (s == i.first)
			return i.second;
	throw std::runtime_error("Invalid argument for option " + name + ". Allowed values are:" + std::accumulate(values.begin(), values.end(), string(), [](const string& s, const pair<string, _t>& v) { return s + ' ' + v.first; }));
}

void Config::set_sens(Sensitivity sens) {
	if (sensitivity != Sensitivity::FAST)
		throw std::runtime_error("Sensitivity switches are mutually exclusive.");
	sensitivity = sens;
}

Config::Config(int argc, const char **argv, bool check_io)
{
	Command_line_parser parser;
	parser.add_command("makedb", "Build DIAMOND database from a FASTA file")
		.add_command("blastp", "Align amino acid query sequences against a protein reference database")
		.add_command("blastx", "Align DNA query sequences against a protein reference database")
		.add_command("view", "View DIAMOND alignment archive (DAA) formatted file")
		.add_command("help", "Produce help message")
		.add_command("version", "Display version information")
		.add_command("getseq", "Retrieve sequences from a DIAMOND database file")
		.add_command("benchmark", "")
		.add_command("random-seqs", "")
		.add_command("compare", "")
		.add_command("sort", "")
		.add_command("roc", "")
		.add_command("dbstat", "")
		.add_command("modelsim", "")
		.add_command("match-file-stat", "")
		.add_command("model-seqs", "")
		.add_command("opt", "")
		.add_command("mask", "")
		.add_command("fastq2fasta", "")
		.add_command("dbinfo", "Print information about a DIAMOND database file")
		.add_command("test-extra", "")
		.add_command("test-io", "")
		.add_command("db-annot-stats", "")
		.add_command("read-sim", "")
		.add_command("info", "")
		.add_command("seed-stat", "")
		.add_command("smith-waterman", "")
		.add_command("cluster", "")
		.add_command("translate", "")
		.add_command("filter-blasttab", "")
		.add_command("show-cbs", "")
		.add_command("simulate-seqs", "")
		.add_command("split", "")
		.add_command("upgma", "")
		.add_command("upgmamc", "")
		.add_command("test", "Run regression tests")
		.add_command("reverse", "")
		.add_command("compute-medoids", "")
		.add_command("mutate", "")
		.add_command("merge-tsv", "");

	string traceback_mode_str;

	Options_group general("General options");
	general.add()
		("threads", 'p', "number of CPU threads", threads_)
		("db", 'd', "database file", database)
		("out", 'o', "output file", output_file)
		("outfmt", 'f', "output format\n\
\t0   = BLAST pairwise\n\
\t5   = BLAST XML\n\
\t6   = BLAST tabular\n\
\t100 = DIAMOND alignment archive (DAA)\n\
\t101 = SAM\n\n\
\tValue 6 may be followed by a space-separated list of these keywords:\n\n\
\tqseqid means Query Seq - id\n\
\tqlen means Query sequence length\n\
\tsseqid means Subject Seq - id\n\
\tsallseqid means All subject Seq - id(s), separated by a ';'\n\
\tslen means Subject sequence length\n\
\tqstart means Start of alignment in query\n\
\tqend means End of alignment in query\n\
\tsstart means Start of alignment in subject\n\
\tsend means End of alignment in subject\n\
\tqseq means Aligned part of query sequence\n\
\tfull_qseq means Query sequence\n\
\tsseq means Aligned part of subject sequence\n\
\tfull_sseq means Subject sequence\n\
\tevalue means Expect value\n\
\tbitscore means Bit score\n\
\tscore means Raw score\n\
\tlength means Alignment length\n\
\tpident means Percentage of identical matches\n\
\tnident means Number of identical matches\n\
\tmismatch means Number of mismatches\n\
\tpositive means Number of positive - scoring matches\n\
\tgapopen means Number of gap openings\n\
\tgaps means Total number of gaps\n\
\tppos means Percentage of positive - scoring matches\n\
\tqframe means Query frame\n\
\tbtop means Blast traceback operations(BTOP)\n\
\tcigar means CIGAR string\n\
\tstaxids means unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)\n\
\tsscinames means unique Subject Scientific Name(s), separated by a ';'\n\
\tsskingdoms means unique Subject Super Kingdom(s), separated by a ';'\n\
\tskingdoms means unique Subject Kingdom(s), separated by a ';'\n\
\tsphylums means unique Subject Phylum(s), separated by a ';'\n\
\tstitle means Subject Title\n\
\tsalltitles means All Subject Title(s), separated by a '<>'\n\
\tqcovhsp means Query Coverage Per HSP\n\
\tscovhsp means Subject Coverage Per HSP\n\
\tqtitle means Query title\n\
\tqqual means Query quality values for the aligned part of the query\n\
\tfull_qqual means Query quality values\n\
\tqstrand means Query strand\n\
\n\tDefault: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", output_format)
("verbose", 'v', "verbose console output", verbose)
("log", 0, "enable debug log", debug_log)
("quiet", 0, "disable console output", quiet)
("header", 0, "Write header lines to blast tabular format.", output_header);

	Options_group makedb("Makedb options");
	makedb.add()
		("in", 0, "input reference file in FASTA format", input_ref_file)
		("taxonmap", 0, "protein accession to taxid mapping file", prot_accession2taxid)
		("taxonnodes", 0, "taxonomy nodes.dmp from NCBI", nodesdmp)
		("taxonnames", 0, "taxonomy names.dmp from NCBI", namesdmp);

	Options_group cluster("");
	cluster.add()
		("cluster-algo", 0, "Clustering algorithm (\"multi-step\", \"mcl\")", cluster_algo)
		("cluster-similarity", 0, "Clustering similarity measure", cluster_similarity)
		("mcl-expansion", 0, "MCL expansion coefficient (default=2)", cluster_mcl_expansion, (uint32_t) 2)
		("mcl-inflation", 0, "MCL inflation coefficient (default=2.0)", cluster_mcl_inflation, 2.0)
		("mcl-sparsity-switch", 0, "MCL switch to sparse matrix computation (default=0.8) ", cluster_mcl_sparsity_switch, 0.8);

	Options_group aligner("Aligner options");
	aligner.add()
		("query", 'q', "input query file", query_file)
		("strand", 0, "query strands to search (both/minus/plus)", query_strands, string("both"))
		("un", 0, "file for unaligned queries", unaligned)
		("al", 0, "file or aligned queries", aligned_file)
		("unfmt", 0, "format of unaligned query file (fasta/fastq)", unfmt, string("fasta"))
		("alfmt", 0, "format of aligned query file (fasta/fastq)", alfmt, string("fasta"))
		("unal", 0, "report unaligned queries (0=no, 1=yes)", report_unaligned, -1)
		("max-target-seqs", 'k', "maximum number of target sequences to report alignments for (default=25)", max_alignments, uint64_t(25))
		("top", 0, "report alignments within this percentage range of top alignment score (overrides --max-target-seqs)", toppercent, 100.0)
		("max-hsps", 0, "maximum number of HSPs per target sequence to report for each query (default=1)", max_hsps, 1u)
		("range-culling", 0, "restrict hit culling to overlapping query ranges", query_range_culling)
		("compress", 0, "compression for output files (0=none, 1=gzip)", compression)
		("evalue", 'e', "maximum e-value to report alignments (default=0.001)", max_evalue, 0.001)
		("min-score", 0, "minimum bit score to report alignments (overrides e-value setting)", min_bit_score)
		("id", 0, "minimum identity% to report an alignment", min_id)
		("query-cover", 0, "minimum query cover% to report an alignment", query_cover)
		("subject-cover", 0, "minimum subject cover% to report an alignment", subject_cover)
		("sensitive", 0, "enable sensitive mode (default: fast)", mode_sensitive)
		("more-sensitive", 0, "enable more sensitive mode (default: fast)", mode_more_sensitive)
		("very-sensitive", 0, "enable very sensitive mode (default: fast)", mode_very_sensitive)
		("ultra-sensitive", 0, "enable ultra sensitive mode (default: fast)", mode_ultra_sensitive)
		("block-size", 'b', "sequence block size in billions of letters (default=2.0)", chunk_size)
		("index-chunks", 'c', "number of chunks for index processing (default=4)", lowmem)
		("tmpdir", 't', "directory for temporary files", tmpdir)
		("parallel-tmpdir", 0, "directory for temporary files used by multiprocessing", parallel_tmpdir)
		("gapopen", 0, "gap open penalty", gap_open, -1)
		("gapextend", 0, "gap extension penalty", gap_extend, -1)
		("frameshift", 'F', "frame shift penalty (default=disabled)", frame_shift)
		("long-reads", 0, "short for --range-culling --top 10 -F 15", long_reads)
		("matrix", 0, "score matrix for protein alignment (default=BLOSUM62)", matrix, string("blosum62"))
		("custom-matrix", 0, "file containing custom scoring matrix", matrix_file)
		("lambda", 0, "lambda parameter for custom matrix", lambda)
		("K", 0, "K parameter for custom matrix", K)
		("comp-based-stats", 0, "enable composition based statistics (0/1=default)", comp_based_stats, 1u)
		("masking", 0, "enable masking of low complexity regions (0/1=default)", masking, 1)
		("query-gencode", 0, "genetic code to use to translate query (see user manual)", query_gencode, 1u)
		("salltitles", 0, "include full subject titles in DAA file", salltitles)
		("sallseqid", 0, "include all subject ids in DAA file", sallseqid)
		("no-self-hits", 0, "suppress reporting of identical self hits", no_self_hits)
		("taxonlist", 0, "restrict search to list of taxon ids (comma-separated)", taxonlist)
		("taxon-exclude", 0, "exclude list of taxon ids (comma-separated)", taxon_exclude);

	Options_group advanced("Advanced options");
	advanced.add()
		("algo", 0, "Seed search algorithm (0=double-indexed/1=query-indexed)", algo, -1)
		("bin", 0, "number of query bins for seed search", query_bins)
		("min-orf", 'l', "ignore translated sequences without an open reading frame of at least this length", run_len)
		("freq-sd", 0, "number of standard deviations for ignoring frequent seeds", freq_sd, 0.0)
		("id2", 0, "minimum number of identities for stage 1 hit", min_identities)
		("xdrop", 'x', "xdrop for ungapped alignment", ungapped_xdrop, 12.3)
		("band", 0, "band for dynamic programming computation", padding)
		("shapes", 's', "number of seed shapes (default=all available)", shapes)
		("shape-mask", 0, "seed shapes", shape_mask)
		("multiprocessing", 0, "enable distributed-memory parallel processing", multiprocessing)
		("mp-init", 0, "initialize multiprocessing run", mp_init)
		("rank-ratio", 0, "include subjects within this ratio of last hit", rank_ratio, -1.0)
		("ext-chunk-size", 0, "chunk size for adaptive ranking (default=400)", ext_chunk_size, (size_t)400)
		("ext", 0, "Extension mode (banded-fast/banded-slow)", ext)
		("culling-overlap", 0, "minimum range overlap with higher scoring hit to delete a hit (default=50%)", inner_culling_overlap, 50.0)
		("taxon-k", 0, "maximum number of targets to report per species", taxon_k, (uint64_t)0)
		("range-cover", 0, "percentage of query range to be covered for range culling (default=50%)", query_range_cover, 50.0)
		("dbsize", 0, "effective database size (in letters)", db_size)
		("no-auto-append", 0, "disable auto appending of DAA and DMND file extensions", no_auto_append)
		("xml-blord-format", 0, "Use gnl|BL_ORD_ID| style format in XML output", xml_blord_format)
		("stop-match-score", 0, "Set the match score of stop codons against each other.", stop_match_score, 1)
		("tantan-minMaskProb", 0, "minimum repeat probability for masking (default=0.9)", tantan_minMaskProb, 0.9)
		("file-buffer-size", 0, "file buffer size in bytes (default=67108864)", file_buffer_size, (size_t)67108864)
		("memory-limit", 'M', "Memory limit for extension stage in GB", memory_limit);

	Options_group view_options("View options");
	view_options.add()
		("daa", 'a', "DIAMOND alignment archive (DAA) file", daa_file)
		("forwardonly", 0, "only show alignments of forward strand", forwardonly);

	Options_group getseq_options("Getseq options");
	getseq_options.add()
		("seq", 0, "Sequence numbers to display.", seq_no);

	double rank_ratio2;
	unsigned window, min_ungapped_score, hit_band, min_hit_score;
	Options_group deprecated_options("");
	deprecated_options.add()
		("window", 'w', "window size for local hit search", window)
		("ungapped-score", 0, "minimum alignment score to continue local extension", min_ungapped_score)
		("hit-band", 0, "band for hit verification", hit_band)
		("hit-score", 0, "minimum score to keep a tentative alignment", min_hit_score)
		("gapped-xdrop", 'X', "xdrop for gapped alignment in bits", gapped_xdrop, 20)
		("rank-ratio2", 0, "include subjects within this ratio of last hit (stage 2)", rank_ratio2, -1.0);

	Options_group hidden_options("");
	hidden_options.add()
		("extend-all", 0, "extend all seed hits", extend_all)
		("local-align", 0, "Local alignment algorithm", local_align_mode, 0u)
		("slow-search", 0, "", slow_search)
		("ht", 0, "", ht_mode)
		("old-freq", 0, "", old_freq)
		("match1", 0, "", match_file1)
		("match2", 0, "", match_file2)
		("max-hits", 'C', "maximum number of hits to consider for one seed", hit_cap)
		("seed-freq", 0, "maximum seed frequency", max_seed_freq, -15.0)
		("space-penalty", 0, "", space_penalty, 0.5)
		("reverse", 0, "", reverse)
		("neighborhood-score", 0, "", neighborhood_score)
		("seed-weight", 'w', "", seed_weight, 7u)
		("idl", 0, "", id_left)
		("idr", 0, "", id_right)
		("idn", 0, "", id_n)
		("bmatch", 0, "", bmatch)
		("bmismatch", 0, "", bmismatch)
		("bcutoff", 0, "", bcutoff)
		("ants", 0, "", n_ants, uint64_t(100))
		("rho", 0, "", rho, 0.99)
		("p_best", 0, "", p_best, 0.05)
		("d_exp", 0, "", d_exp, 1.0)
		("d_new", 0, "", d_new, 1.0)
		("score-estimate-factor", 0, "", score_estimate_factor, 0.0)
		("diag-min-estimate", 0, "", diag_min_estimate, 17)
		("qfilt", 0, "", qfilt)
		("sfilt", 0, "", sfilt)
		("path-cutoff", 0, "", path_cutoff, 0.92)
		("sw", 0, "", use_smith_waterman)
		("superblock", 0, "", superblock, 128)
		("max-cells", 0, "", max_cells, 10000000u)
		("load-balancing", 0, "", load_balancing, (unsigned)Config::query_parallel)
		("br", 0, "", benchmark_ranking)
		("log-query", 0, "", log_query)
		("log-subject", 0, "", log_subject)
		("palign", 0, "", threads_align)
		("score-ratio", 0, "", score_ratio, 0.9)
		("fetch-size", 0, "trace point fetch size", fetch_size, 4096u)
		("target-fetch-size", 0, "number of target sequences to fetch for seed extension", target_fetch_size, 4u)
		("small-query", 0, "", small_query)
		("hashed-seeds", 0, "", hashed_seeds)
		("rank-factor", 0, "", rank_factor, -1.0)
		("transcript-len-estimate", 0, "", transcript_len_estimate, 1.0)
		("family-counts", 0, "", family_counts_file)
		("radix-cluster-buffered", 0, "", radix_cluster_buffered)
		("join-split-size", 0, "", join_split_size, 100000u)
		("join-split-key-len", 0, "", join_split_key_len, 17u)
		("radix-bits", 0, "", radix_bits, 8u)
		("join-ht-factor", 0, "", join_ht_factor, 1.3)
		("sort-join", 0, "", sort_join)
		("simple-freq", 0, "", simple_freq)
		("freq-treshold", 0, "", freq_treshold)
		("filter-locus", 0, "", filter_locus)
		("use-dataset-field", 0, "", use_dataset_field)
		("store-query-quality", 0, "", store_query_quality)
		("swipe-chunk-size", 0, "", swipe_chunk_size, 256u)
		("query-parallel-limit", 0, "", query_parallel_limit, 3000000u)
		("hard-masked", 0, "", hardmasked)
		("cbs-window", 0, "", cbs_window, 40)
		("no-unlink", 0, "", no_unlink)
		("no-dict", 0, "", no_dict)
		("swipe", 0, "", swipe_all)
		("upgma-edge-limit", 0, "", upgma_edge_limit, (uint64_t)10000000)
		("tree", 0, "", tree_file)
		("upgma-dist", 0, "", upgma_dist)
		("upgma-input", 0, "", upgma_input)
		("log-extend", 0, "", log_extend)
		("chaining-maxgap", 0, "", chaining_maxgap, 2000)
		("tantan-maxRepeatOffset", 0, "maximum tandem repeat period to consider (50)", tantan_maxRepeatOffset, 15)
		("tantan-ungapped", 0, "use tantan masking in ungapped mode", tantan_ungapped)
		("family-map", 0, "", family_map)
		("family-map-query", 0, "", family_map_query)
		("chaining-range-cover", 0, "", chaining_range_cover, (size_t)8)
		("index-mode", 0, "index mode (0=4x12, 1=16x9)", index_mode)
		("no-swipe-realign", 0, "", no_swipe_realign)
		("cut-bar", 0, "", cut_bar)
		("bootstrap", 0, "", bootstrap)
		("chaining-maxnodes", 0, "", chaining_maxnodes)
		("cutoff-score-8bit", 0, "", cutoff_score_8bit, 240)
		("min-band-overlap", 0, "", min_band_overlap, 0.2)
		("min-realign-overhang", 0, "", min_realign_overhang, 30)
		("ungapped-window", 0, "", ungapped_window, 48)
		("gapped-filter-diag-score", 0, "", gapped_filter_diag_score, 20)
		("gapped-filter-evalue", 0, "", gapped_filter_evalue, -1.0)
		("gapped-filter-window", 0, "", gapped_filter_window, 200)
		("output-hits", 0, "", output_hits)
		("ungapped-evalue", 0, "", ungapped_evalue, -1.0)
		("no-logfile", 0, "", no_logfile)
		("no-heartbeat", 0, "", no_heartbeat)
		("band-bin", 0, "", band_bin, 24)
		("col-bin", 0, "", col_bin, 400)
		("self", 0, "", self)
		("trace-pt-fetch-size", 0, "", trace_pt_fetch_size, (size_t)10e9)
		("tile-size", 0, "", tile_size, (uint32_t)1024)
		("short-query-ungapped-bitscore", 0, "", short_query_ungapped_bitscore, 25.0)
		("short-query-max-len", 0, "", short_query_max_len, 60)
		("gapped-filter-evalue1", 0, "", gapped_filter_evalue1, 1.0e+04)
		("ext-yield", 0, "", ext_min_yield)
		("full-sw-len", 0, "", full_sw_len)
		("relaxed-evalue-factor", 0, "", relaxed_evalue_factor, 1.0)
		("type", 0, "", type)
		("raw", 0, "", raw)
		("chaining-len-cap", 0, "", chaining_len_cap, 2.0)
		("chaining-min-nodes", 0, "", chaining_min_nodes, (size_t)200)
		("fast-tsv", 0, "", fast_tsv)
		("target-parallel-verbosity", 0, "", target_parallel_verbosity, UINT_MAX)
		("ext-targets", 0, "", global_ranking_targets)
		("traceback-mode", 0, "", traceback_mode_str)
		("mid-sensitive", 0, "", mode_mid_sensitive);
	
	parser.add(general).add(makedb).add(cluster).add(aligner).add(advanced).add(view_options).add(getseq_options).add(hidden_options).add(deprecated_options);
	parser.store(argc, argv, command);

	traceback_mode = set_string_option<TracebackMode>(traceback_mode_str, "--traceback-mode",
		{ {"score", TracebackMode::SCORE_ONLY },
		{"stat", TracebackMode::STAT},
		{"vector", TracebackMode::VECTOR},
		{"buffer", TracebackMode::SCORE_BUFFER} });

	if (toppercent != 100.0 && max_alignments != 25)
		throw std::runtime_error("--top and --max-target-seqs are mutually exclusive.");

	if (long_reads) {
		query_range_culling = true;
		if (toppercent == 100.0)
			toppercent = 10.0;
		if (frame_shift == 0)
			frame_shift = 15;
	}

	if (check_io) {
		switch (command) {
		case Config::makedb:
			if (database == "")
				throw std::runtime_error("Missing parameter: database file (--db/-d)");
			if (chunk_size != 0.0)
				throw std::runtime_error("Invalid option: --block-size/-b. Block size is set for the alignment commands.");
			break;
		case Config::blastp:
		case Config::blastx:
			if (database == "")
				throw std::runtime_error("Missing parameter: database file (--db/-d)");
			if (daa_file.length() > 0) {
				if (output_file.length() > 0)
					throw std::runtime_error("Options --daa and --out cannot be used together.");
				if (output_format.size() > 0 && output_format[0] != "daa")
					throw std::runtime_error("Invalid parameter: --daa/-a. Output file is specified with the --out/-o parameter.");
				output_file = daa_file;
			}
			if (daa_file.length() > 0 || (output_format.size() > 0 && (output_format[0] == "daa" || output_format[0] == "100"))) {
				if (compression != 0)
					throw std::runtime_error("Compression is not supported for DAA format.");
				if (!no_auto_append)
					auto_append_extension(output_file, ".daa");
			}
			break;
		case Config::view:
			if (daa_file == "")
				throw std::runtime_error("Missing parameter: DAA file (--daa/-a)");
		default:
			;
		}

		switch (command) {
		case Config::dbinfo:
			if (database == "")
				throw std::runtime_error("Missing parameter: database file (--db/-d)");
		}
	}

	if (hit_cap != 0)
		throw std::runtime_error("Deprecated parameter: --max-hits/-C.");

	if (debug_log)
		verbosity = 3;
	else if (quiet)
		verbosity = 0;
	else if (verbose)
		verbosity = 2;
	else if (((command == Config::view || command == blastx || command == blastp) && output_file == "")
		|| command == Config::version || command == getseq || command == fastq2fasta || command == regression_test)
		verbosity = 0;
	else
		verbosity = 1;

	switch (verbosity) {
	case 0:
		message_stream = Message_stream(false);
		break;
	case 3:
		log_stream = Message_stream(true, !config.no_logfile);
		verbose_stream = Message_stream(true, !config.no_logfile);
		message_stream = Message_stream(true, !config.no_logfile);
		break;
	case 2:
		verbose_stream = Message_stream();
	default:
		;
	}

	invocation = join(" ", vector<string>(&argv[0], &argv[argc]));
	log_stream << invocation << endl;

	if (!no_auto_append) {
		if (command == Config::makedb)
			auto_append_extension(database, ".dmnd");
		else
			auto_append_extension_if_exists(database, ".dmnd");
		if (command == Config::view)
			auto_append_extension(daa_file, ".daa");
		if (compression == 1)
			auto_append_extension(output_file, ".gz");
	}

	if (verbosity >= 1 || command == regression_test) {
		ostream &header_out = command == Config::help ? cout : cerr;
		header_out << Const::program_name << " v" << Const::version_string << "." << (unsigned)Const::build_version << " (C) Max Planck Society for the Advancement of Science" << endl;
		header_out << "Documentation, support and updates available at http://www.diamondsearch.org" << endl << endl;
	}
	log_stream << Const::program_name << " v" << Const::version_string << "." << (unsigned)Const::build_version << endl;
#ifndef NDEBUG
	verbose_stream << "Assertions enabled." << endl;
#endif
	set_option(threads_, std::thread::hardware_concurrency());

	switch (command) {
	case Config::makedb:
	case Config::blastp:
	case Config::blastx:
	case Config::view:
	case Config::cluster:
	case Config::regression_test:
	case Config::compute_medoids:
		message_stream << "#CPU threads: " << threads_ << endl;
	default:
		;
	}

	switch (command) {
	case Config::blastp:
	case Config::blastx:
	case Config::benchmark:
	case Config::model_sim:
	case Config::opt:
	case Config::mask:
	case Config::makedb:
	case Config::cluster:
	case Config::regression_test:
	case Config::compute_medoids:
		if (frame_shift != 0 && command == Config::blastp)
			throw std::runtime_error("Frameshift alignments are only supported for translated searches.");
		if (query_range_culling && frame_shift == 0)
			throw std::runtime_error("Query range culling is only supported in frameshift alignment mode (option -F).");
		if (matrix_file == "")
			score_matrix = Score_matrix(to_upper_case(matrix), gap_open, gap_extend, frame_shift, stop_match_score);
		else {
			if (lambda == 0 || K == 0)
				throw std::runtime_error("Custom scoring matrices require setting the --lambda and --K options.");
			if (gap_open == -1 || gap_extend == -1)
				throw std::runtime_error("Custom scoring matrices require setting the --gapopen and --gapextend options.");
			score_matrix = Score_matrix(matrix_file, lambda, K, gap_open, gap_extend);
		}
		if(command == Config::cluster && !Workflow::Cluster::ClusterRegistry::has(cluster_algo)){
			ostream &header_out = command == Config::help ? cout : cerr;
			header_out << "Unkown clustering algorithm: " << cluster_algo << endl;
			header_out << "Available options are: " << endl;
			for(string c_algo : Workflow::Cluster::ClusterRegistry::getKeys()){
				header_out << "\t" << c_algo << "\t"<< Workflow::Cluster::ClusterRegistry::get(c_algo)->get_description() << endl;
			}
			throw std::runtime_error("Clustering algorithm not found.");
		}
		message_stream << "Scoring parameters: " << score_matrix << endl;
		if (masking == 1)
			Masking::instance = unique_ptr<Masking>(new Masking(score_matrix));
	}

	if (command == Config::blastp || command == Config::blastx || command == Config::benchmark || command == Config::model_sim || command == Config::opt
		|| command == Config::mask || command == Config::cluster || command == Config::compute_medoids || command == Config::regression_test) {
		if (tmpdir == "")
			tmpdir = extract_dir(output_file);

		init_cbs();
		raw_ungapped_xdrop = score_matrix.rawscore(ungapped_xdrop);
		verbose_stream << "CPU features detected: " << SIMD::features() << endl;
	}

	sensitivity = Sensitivity::FAST;
	if (mode_mid_sensitive) set_sens(Sensitivity::MID_SENSITIVE);
	if (mode_sensitive) set_sens(Sensitivity::SENSITIVE);
	if (mode_more_sensitive) set_sens(Sensitivity::MORE_SENSITIVE);
	if (mode_very_sensitive) set_sens(Sensitivity::VERY_SENSITIVE);
	if (mode_ultra_sensitive) set_sens(Sensitivity::ULTRA_SENSITIVE);

	if (ext != "banded-fast" && ext != "banded-slow" && ext != "")
		throw std::runtime_error("Possible values for --ext are: banded-fast, banded-slow");

	Translator::init(query_gencode);

	if (command == blastx)
		input_value_traits = nucleotide_traits;

	if (command == help)
		parser.print_help();

	if (query_strands != "both" && query_strands != "minus" && query_strands != "plus")
		throw std::runtime_error("Invalid value for parameter --strand");

	if (unfmt == "fastq" || alfmt == "fastq")
		store_query_quality = true;
	if (!aligned_file.empty())
		log_stream << "Aligned file format: " << alfmt << endl;

	/*log_stream << "sizeof(hit)=" << sizeof(hit) << " sizeof(packed_uint40_t)=" << sizeof(packed_uint40_t)
		<< " sizeof(sorted_list::entry)=" << sizeof(sorted_list::entry) << endl;*/

	if (swipe_all) {
		algo = double_indexed;
	}

	use_lazy_dict = false;

	if (query_range_culling && taxon_k != 0)
		throw std::runtime_error("--taxon-k is not supported for --range-culling mode.");

	if (parallel_tmpdir == "") {
		parallel_tmpdir = tmpdir;
	} else {
#ifndef WIN32
		if (multiprocessing) {
			// char * env_str = std::getenv("SLURM_JOBID");
			// if (env_str) {
			// 	parallel_tmpdir = join_path(parallel_tmpdir, "diamond_job_"+string(env_str));
			// }
			errno = 0;
			int s = mkdir(parallel_tmpdir.c_str(), 00770);
			if (s != 0) {
				if (errno == EEXIST) {
					// directory did already exist
				} else {
					throw(std::runtime_error("could not create parallel temporary directory " + parallel_tmpdir));
				}
			}
		}
#endif
	}

	log_stream << "MAX_SHAPE_LEN=" << MAX_SHAPE_LEN;
#ifdef SEQ_MASK
	log_stream << " SEQ_MASK";
#endif
#ifdef STRICT_BAND
	log_stream << " STRICT_BAND";
#endif
	log_stream << endl;
}
