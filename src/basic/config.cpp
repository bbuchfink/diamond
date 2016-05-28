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

Config config;

Config::Config(int argc, const char **argv)
{
	Command_line_parser parser;
	parser.add_command("makedb", "Build DIAMOND database from a FASTA file")
	      .add_command("blastp", "Align amino acid query sequences against a protein reference database")
		  .add_command("blastx", "Align DNA query sequences against a protein reference database")
		  .add_command("view", "View DIAMOND alignment archive (DAA) formatted file")
		  .add_command("help", "Produce help message")
		  .add_command("version", "Display version information");

	Options_group general ("General options");
	general.add()
		("threads", 'p', "number of CPU threads", threads_)
		("db", 'd', "database file", database)
		("daa", 'a', "DIAMOND alignment archive (DAA) file", daa_file)
		("verbose", 'v', "verbose console output", verbose)
		("log", 0, "enable debug log", debug_log)
		("quiet", 0, "disable console output", quiet);

	Options_group makedb("Makedb options");
	makedb.add()
		("in", 0, "input reference file in FASTA format", input_ref_file)
		("block-size", 'b', "sequence block size in billions of letters (default=2)", chunk_size)
#ifdef EXTRA
		("dbtype", po::value<string>(&program_options::db_type), "database type (nucl/prot)")
#endif
		;

	Options_group aligner("Aligner options");
	aligner.add()
		("query",'q', "input query file", query_file)
		("max-target-seqs",'k', "maximum number of target sequences to report alignments for", max_alignments, 25lu)
		("top", 0, "report alignments within this percentage range of top alignment score (overrides --max-target-seqs)", toppercent, 100.0)
		("compress", 0, "compression for output files (0=none, 1=gzip)", compression)
		("evalue",'e', "maximum e-value to report alignments", max_evalue, 0.001)
		("min-score", 0, "minimum bit score to report alignments (overrides e-value setting)", min_bit_score)
		("id", 0, "minimum identity% to report an alignment", min_id)
		("query-cover", 0, "minimum query cover% to report an alignment", query_cover)
		("sensitive", 0, "enable sensitive mode (default: fast)", mode_sensitive)
		("index-chunks",'c', "number of chunks for index processing", lowmem, 4u)
		("tmpdir",'t', "directory for temporary files", tmpdir)
		("gapopen", 0, "gap open penalty (default=11 for protein)", gap_open, -1)
		("gapextend", 0, "gap extension penalty (default=1 for protein)", gap_extend, -1)
#ifdef EXTRA
		("reward", po::value<int>(&program_options::reward)->default_value(2), "match reward score (blastn only)")
		("penalty", po::value<int>(&program_options::penalty)->default_value(-3), "mismatch penalty score (blastn only)")
#endif
		("matrix", 0, "score matrix for protein alignment", matrix, string("blosum62"))
		("seg", 0, "enable SEG masking of queries (yes/no)", seg)
		("salltitles", 0, "print full subject titles in output files", salltitles);

	Options_group advanced("Advanced options");
	advanced.add()
		("seed-freq", 0, "maximum seed frequency", max_seed_freq, -15.0)
		("run-len",'l', "mask runs between stop codons shorter than this length", run_len)
		("max-hits",'C', "maximum number of hits to consider for one seed", hit_cap)
		("id2", 0, "minimum number of identities for stage 1 hit", min_identities)
		("window", 'w', "window size for local hit search", window)
		("xdrop", 'x', "xdrop for ungapped alignment", xdrop, 20)
		("gapped-xdrop",'X', "xdrop for gapped alignment in bits", gapped_xdrop, 20)
		("ungapped-score", 0, "minimum raw alignment score to continue local extension", min_ungapped_raw_score)
		("hit-band", 0, "band for hit verification", hit_band)
		("hit-score",0, "minimum score to keep a tentative alignment", min_hit_score)
		("band", 0, "band for dynamic programming computation", padding)
		("local-align", 0, "Local alignment algorithm", local_align_mode, 0u)
		("shapes", 's', "number of seed shapes (0 = all available)", shapes)
		("index-mode", 0, "index mode (0=4x12, 1=16x9)", index_mode)
		("fetch-size", 0, "trace point fetch size", fetch_size, 4096u)
		("single-domain", 0, "Discard secondary domains within one target sequence", single_domain)
		("dbsize", 0, "effective database size (in letters)", db_size)
		("no-auto-append", 0, "disable auto appending of DAA and DMND file extensions", no_auto_append);
	
	Options_group view_options("View options");
	view_options.add()
		("out", 'o', "output file", output_file)
		("outfmt",'f', "output format (tab/sam)", output_format, string("tab"))
		("forwardonly", 0, "only show alignments of forward strand", forwardonly);

#ifdef EXTRA
	("match1", po::value<string>(&program_options::match_file1))
		("match2", po::value<string>(&program_options::match_file2))
		("tab", "tabular format")
#endif

	parser.add(general).add(makedb).add(aligner).add(advanced).add(view_options);
	parser.store(argc, argv, command);

	if (debug_log)
		verbosity = 3;
	else if (quiet)
		verbosity = 0;
	else if (verbose)
		verbosity = 2;
	else if ((command == Config::view && output_file == "")
		|| command == Config::version)
		verbosity = 0;
	else
		verbosity = 1;

	switch (verbosity) {
	case 0:
		message_stream = Message_stream(false);
		break;	
	case 3:
		log_stream = Message_stream();
	case 2:
		verbose_stream = Message_stream();
	default:
		;
	}

	switch (command) {
	case Config::makedb:
		if (input_ref_file == "")
			throw std::runtime_error("Missing parameter: input file (--in)");
		if (database == "")
			throw std::runtime_error("Missing parameter: database file (--db/-d)");
		set_option(chunk_size, 2.0);
		break;
	case Config::blastp:
	case Config::blastx:
		if (query_file == "")
			throw std::runtime_error("Missing parameter: query file (--query/-q)");
		if (database == "")
			throw std::runtime_error("Missing parameter: database file (--db/-d)");
		if (chunk_size != 0)
			std::cerr << "Warning: --block-size option should be set for the makedb command." << endl;
	case Config::view:
		if (daa_file == "")
			throw std::runtime_error("Missing parameter: DAA file (--daa/-a)");
	default:
		;
	}
	
	if (mode_sensitive) {
		set_option(index_mode, 1u);
		//lowmem = std::max(lowmem, 4u);
	}
	else {
		set_option(index_mode, 0u);
	}

	if (!no_auto_append) {
		auto_append_extension(database, ".dmnd");
		auto_append_extension(daa_file, ".daa");
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

	if (command == Config::blastp || command == Config::blastx) {
		if (tmpdir == "")
			tmpdir = extract_dir(daa_file);
		if (gap_open == -1)
			gap_open = 11;
		if (gap_extend == -1)
			gap_extend = 1;
		score_matrix = Score_matrix(matrix, gap_open, gap_extend, reward, penalty);
		message_stream << "Scoring parameters: " << score_matrix << endl;
		if (seg == "" && command == blastx)
			seg = "yes";
		verbose_stream << "SEG masking = " << (seg == "yes") << endl;
		have_ssse3 = check_SSSE3();
		if (have_ssse3)
			verbose_stream << "SSSE3 enabled." << endl;
		verbose_stream << "Reduction: " << Reduction::reduction << endl;
		::shapes = shape_config(index_mode, shapes);
		verbose_stream << "Shape configuration: " << ::shapes << endl;
		message_stream << "#Target sequences to report alignments for: ";
		if (max_alignments == 0) {
			max_alignments = std::numeric_limits<unsigned long>::max();
			message_stream << "unlimited" << endl;
		} else
			message_stream << max_alignments << endl;
	}

	if (command == blastx)
		input_value_traits = nucleotide_traits;
	if (command == help)
		parser.print_help();

	log_stream << "sizeof(hit)=" << sizeof(hit) << " sizeof(packed_uint40_t)=" << sizeof(packed_uint40_t)
		<< " sizeof(sorted_list::entry)=" << sizeof(sorted_list::entry) << endl;
}

void Config::set_chunk_size(double x)
{
	set_option(chunk_size, x);
}
