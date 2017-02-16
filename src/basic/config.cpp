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

#include "../util/command_line_parser.h"
#include "config.h"
#include "../util/util.h"
#include "../util/log_stream.h"
#include "../util/tinythread.h"
#include "../basic/value.h"
#include "score_matrix.h"
#include "../util/system.h"
#include "reduction.h"
#include "shape_config.h"
#include "../util/temp_file.h"
#include "../basic/match.h"
#include "../data/sorted_list.h"
#include "../basic/translate.h"
#include "../dp/dp.h"

Config config;

Config::Config(int argc, const char **argv)
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
		.add_command("opt", "");

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
\tsseq means Aligned part of subject sequence\n\
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
\tstitle means Subject Title\n\
\tsalltitles means All Subject Title(s), separated by a '<>'\n\
\tqcovhsp means Query Coverage Per HSP\n\
\tqtitle means Query title\n\n\
\tDefault: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", output_format)
("verbose", 'v', "verbose console output", verbose)
("log", 0, "enable debug log", debug_log)
("quiet", 0, "disable console output", quiet);

	Options_group makedb("Makedb options");
	makedb.add()
		("in", 0, "input reference file in FASTA format", input_ref_file);

	Options_group aligner("Aligner options");
	aligner.add()
		("query", 'q', "input query file", query_file)
		("un", 0, "file for unaligned queries", unaligned)
		("unal", 0, "report unaligned queries (0=no, 1=yes)", report_unaligned, -1)
		("max-target-seqs", 'k', "maximum number of target sequences to report alignments for", max_alignments, uint64_t(25))
		("top", 0, "report alignments within this percentage range of top alignment score (overrides --max-target-seqs)", toppercent, 100.0)
		("compress", 0, "compression for output files (0=none, 1=gzip)", compression)
		("evalue", 'e', "maximum e-value to report alignments (default=0.001)", max_evalue, 0.001)
		("min-score", 0, "minimum bit score to report alignments (overrides e-value setting)", min_bit_score)
		("id", 0, "minimum identity% to report an alignment", min_id)
		("query-cover", 0, "minimum query cover% to report an alignment", query_cover)
		("subject-cover", 0, "minimum subject cover% to report an alignment", subject_cover)
		("sensitive", 0, "enable sensitive mode (default: fast)", mode_sensitive)
		("more-sensitive", 0, "enable more sensitive mode (default: fast)", mode_more_sensitive)
		("block-size", 'b', "sequence block size in billions of letters (default=2.0)", chunk_size)
		("index-chunks", 'c', "number of chunks for index processing", lowmem)
		("tmpdir", 't', "directory for temporary files", tmpdir)
		("gapopen", 0, "gap open penalty (default=11 for protein)", gap_open, -1)
		("gapextend", 0, "gap extension penalty (default=1 for protein)", gap_extend, -1)
		("matrix", 0, "score matrix for protein alignment (default=BLOSUM62)", matrix, string("blosum62"))
		("custom-matrix", 0, "file containing custom scoring matrix", matrix_file)
		("lambda", 0, "lambda parameter for custom matrix", lambda)
		("K", 0, "K parameter for custom matrix", K)
		("comp-based-stats", 0, "enable composition based statistics (0/1=default)", comp_based_stats, 1u)
		("seg", 0, "enable SEG masking of queries (yes/no)", seg)
		("query-gencode", 0, "genetic code to use to translate query (see user manual)", query_gencode, 1u)
		("salltitles", 0, "print full subject titles in output files", salltitles)
		("no-self-hits", 0, "suppress reporting of identical self hits", no_self_hits);

	Options_group advanced("Advanced options");
	advanced.add()
		("bin", 0, "number of query bins for seed search", query_bins, 16u)
		("min-orf", 'l', "ignore translated sequences without an open reading frame of at least this length", run_len)
		("freq-sd", 0, "number of standard deviations for ignoring frequent seeds", freq_sd, 0.0)
		("id2", 0, "minimum number of identities for stage 1 hit", min_identities)
		("window", 'w', "window size for local hit search", window)
		("xdrop", 'x', "xdrop for ungapped alignment", ungapped_xdrop, 12.3)
		("ungapped-score", 0, "minimum alignment score to continue local extension", min_ungapped_score)
		("hit-band", 0, "band for hit verification", hit_band)
		("hit-score", 0, "minimum score to keep a tentative alignment", min_hit_score)
		("gapped-xdrop", 'X', "xdrop for gapped alignment in bits", gapped_xdrop, 20)
		("band", 0, "band for dynamic programming computation", padding)
		("shapes", 's', "number of seed shapes (0 = all available)", shapes)
		("shape-mask", 0, "seed shapes", shape_mask)
		("index-mode", 0, "index mode (0=4x12, 1=16x9)", index_mode)
		("fetch-size", 0, "trace point fetch size", fetch_size, 4096u)
		("rank-factor", 0, "include subjects within this range of max-target-seqs", rank_factor, 2.0)
		("rank-ratio", 0, "include subjects within this ratio of last hit", rank_ratio, 0.6)
		("max-hsps", 0, "maximum number of HSPs per subject sequence to save for each query", max_hsps, 1u)
		("dbsize", 0, "effective database size (in letters)", db_size)
		("no-auto-append", 0, "disable auto appending of DAA and DMND file extensions", no_auto_append)
		("target-fetch-size", 0, "number of target sequences to fetch for seed extension", target_fetch_size, 4u);

	Options_group view_options("View options");
	view_options.add()
		("daa", 'a', "DIAMOND alignment archive (DAA) file", daa_file)
		("forwardonly", 0, "only show alignments of forward strand", forwardonly);

	Options_group getseq_options("Getseq options");
	getseq_options.add()
		("seq", 0, "Sequence numbers to display.", seq_no);

	Options_group hidden_options("");
	hidden_options.add()
		("extend-all", 0, "extend all seed hits", extend_all)
		("local-align", 0, "Local alignment algorithm", local_align_mode, 0u)
		("slow-search", 0, "", slow_search)
		("ht", 0, "", ht_mode)
		("old-freq", 0, "", old_freq)
		("qp", 0, "", query_parallel)
		("match1", 0, "", match_file1)
		("match2", 0, "", match_file2)
		("max-hits", 'C', "maximum number of hits to consider for one seed", hit_cap)
		("seed-freq", 0, "maximum seed frequency", max_seed_freq, -15.0)
		("space-penalty", 0, "", space_penalty, 4.9)
		("min-diag-score", 0, "", min_diag_score, 10.5)
		("reverse", 0, "", reverse)
		("neighborhood-score", 0, "", neighborhood_score)
		("algo", 0, "", algo, 0u)
		("seed-weight", 'w', "", seed_weight, 7u)
		("very-sensitive", 0, "", mode_very_sensitive)
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
		("greedy", 0, "", greedy)
		("qfilt", 0, "", qfilt);
		
	parser.add(general).add(makedb).add(aligner).add(advanced).add(view_options).add(getseq_options).add(hidden_options);
	parser.store(argc, argv, command);

	switch (command) {
	case Config::makedb:
		if (input_ref_file == "")
			throw std::runtime_error("Missing parameter: input file (--in)");
		if (database == "")
			throw std::runtime_error("Missing parameter: database file (--db/-d)");
		if (chunk_size != 0.0)
			throw std::runtime_error("Invalid option: --block-size/-b. Block size is set for the alignment commands.");
		break;
	case Config::blastp:
	case Config::blastx:
		if (query_file == "")
			throw std::runtime_error("Missing parameter: query file (--query/-q)");
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

	if (hit_cap != 0)
		throw std::runtime_error("Deprecated parameter: --max-hits/-C.");

	if (debug_log)
		verbosity = 3;
	else if (quiet)
		verbosity = 0;
	else if (verbose)
		verbosity = 2;
	else if (((command == Config::view || command == blastx || command == blastp) && output_file == "")
		|| command == Config::version || command == getseq)
		verbosity = 0;
	else
		verbosity = 1;

	switch (verbosity) {
	case 0:
		message_stream = Message_stream(false);
		break;
	case 3:
		log_stream = Message_stream(true, true);
		verbose_stream = Message_stream(true, true);
		message_stream = Message_stream(true, true);
		break;
	case 2:
		verbose_stream = Message_stream();
	default:
		;
	}

	for (int i = 0; i < argc; ++i)
		log_stream << argv[i] << ' ';
	log_stream << endl;

	if (!no_auto_append) {
		auto_append_extension(database, ".dmnd");
		if (command == Config::view)
			auto_append_extension(daa_file, ".daa");
		if (compression == 1)
			auto_append_extension(output_file, ".gz");
	}

	message_stream << Const::program_name << " v" << Const::version_string << "." << (unsigned)Const::build_version << " | by Benjamin Buchfink <buchfink@gmail.com>" << endl;
	message_stream << "Check http://github.com/bbuchfink/diamond for updates." << endl << endl;
#ifndef NDEBUG
	verbose_stream << "Assertions enabled." << endl;
#endif
	set_option(threads_, tthread::thread::hardware_concurrency());

	switch (command) {
	case Config::makedb:
	case Config::blastp:
	case Config::blastx:
	case Config::view:
		message_stream << "#CPU threads: " << threads_ << endl;
	default:
		;
	}

	if (command == Config::blastp || command == Config::blastx || command == Config::benchmark || command == Config::model_sim || command == Config::opt) {
		if (tmpdir == "")
			tmpdir = extract_dir(output_file);
		if (gap_open == -1)
			gap_open = 11;
		if (gap_extend == -1)
			gap_extend = 1;
		if (matrix_file == "")
			score_matrix = Score_matrix(to_upper_case(matrix), gap_open, gap_extend, reward, penalty);
		else {
			if (lambda == 0 || K == 0)
				throw std::runtime_error("Custom scoring matrices require setting the --lambda and --K options.");
			score_matrix = Score_matrix(matrix_file, lambda, K, gap_open, gap_extend);
		}
		message_stream << "Scoring parameters: " << score_matrix << endl;
		raw_ungapped_xdrop = score_matrix.rawscore(ungapped_xdrop);
		min_diag_raw_score = score_matrix.rawscore(min_diag_score);
		raw_space_penalty = score_matrix.rawscore(space_penalty, double());
		log_stream << "Min_diag_score=" << min_diag_raw_score << " space_penalty=" << raw_space_penalty << endl;
		init_cbs();

		if (seg == "" && command == blastx)
			seg = "yes";
		verbose_stream << "SEG masking = " << (seg == "yes") << endl;
		have_ssse3 = check_SSSE3();
		if (have_ssse3)
			verbose_stream << "SSSE3 enabled." << endl;
		verbose_stream << "Reduction: " << Reduction::reduction << endl;

		if (mode_more_sensitive) {
			set_option(index_mode, 1u);
			set_option(freq_sd, 200.0);
		}
		else if (mode_sensitive) {
			set_option(index_mode, 1u);
			set_option(freq_sd, 10.0);
		}
		else {
			set_option(index_mode, 0u);
			set_option(freq_sd, 50.0);
		}

		verbose_stream << "Seed frequency SD: " << freq_sd << endl;
		::shapes = shape_config(index_mode, shapes, shape_mask);
		verbose_stream << "Shape configuration: " << ::shapes << endl;
		seed_anchor = std::min(::shapes[0].length_ - 1, 8u);

		message_stream << "#Target sequences to report alignments for: ";
		if (max_alignments == 0) {
			max_alignments = std::numeric_limits<uint64_t>::max();
			message_stream << "unlimited" << endl;
		}
		else
			message_stream << max_alignments << endl;
	}

	Translator::init(query_gencode);

	if (command == blastx)
		input_value_traits = nucleotide_traits;

	if (command == help)
		parser.print_help();

	/*log_stream << "sizeof(hit)=" << sizeof(hit) << " sizeof(packed_uint40_t)=" << sizeof(packed_uint40_t)
		<< " sizeof(sorted_list::entry)=" << sizeof(sorted_list::entry) << endl;*/
}
