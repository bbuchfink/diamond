 /****
DIAMOND protein aligner
Copyright (C) 2013-2022 Max Planck Society for the Advancement of Science e.V.
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
#include <iostream>
#include <cstdlib>
#include <sys/stat.h>
#include <exception>
#include <iomanip>
#include <numeric>
#include <thread>
#include "../util/command_line_parser.h"
#include "config.h"
#include "../util/util.h"
#include "../util/log_stream.h"
#include "../basic/value.h"
#include "../stats/score_matrix.h"
#include "../util/system.h"
#include "reduction.h"
#include "shape_config.h"
#include "../util/io/temp_file.h"
#include "../basic/match.h"
#include "../cluster/cluster_registry.h"
#include "../util/sequence/translate.h"
#include "../dp/dp.h"
#include "../masking/masking.h"
#include "../util/system/system.h"
#include "../util/simd.h"
#include "../util/parallel/multiprocessing.h"
#include "../search/search.h"

using std::thread;
using std::stringstream;
using std::endl;
using std::ostream;
using std::unique_ptr;
using std::string;
using std::pair;
using std::cout;
using std::cerr;

const EMap<Sensitivity> EnumTraits<Sensitivity>::to_string = {
	{ Sensitivity::FASTER, "faster" },
	{ Sensitivity::FAST, "fast" },
	{ Sensitivity::DEFAULT, "default" },
	{ Sensitivity::MID_SENSITIVE, "mid-sensitive" },
	{ Sensitivity::SENSITIVE, "sensitive" },
	{ Sensitivity::MORE_SENSITIVE, "more-sensitive" },
	{ Sensitivity::VERY_SENSITIVE, "very-sensitive" },
	{ Sensitivity::ULTRA_SENSITIVE, "ultra-sensitive" }
};

const SEMap<Sensitivity> EnumTraits<Sensitivity>::from_string = {
	{ "faster", Sensitivity::FASTER },
	{ "fast", Sensitivity::FAST },
	{ "default", Sensitivity::DEFAULT },
	{ "mid-sensitive", Sensitivity::MID_SENSITIVE },
	{ "sensitive", Sensitivity::SENSITIVE },
	{ "more-sensitive", Sensitivity::MORE_SENSITIVE },
	{ "very-sensitive", Sensitivity::VERY_SENSITIVE },
	{ "ultra-sensitive", Sensitivity::ULTRA_SENSITIVE }
};

const SEMap<GraphAlgo> EnumTraits<GraphAlgo>::from_string = {
	{ "gvc", GraphAlgo::GREEDY_VERTEX_COVER },
	{ "len", GraphAlgo::LEN_SORTED }
};

const EMap<Config::Algo> EnumTraits<Config::Algo>::to_string = { { Config::Algo::DOUBLE_INDEXED, "Double-indexed" }, { Config::Algo::QUERY_INDEXED, "Query-indexed"}, {Config::Algo::CTG_SEED, "Query-indexed with contiguous seed"} };
const SEMap<Config::Algo> EnumTraits<Config::Algo>::from_string = { {"", Config::Algo::AUTO}, { "0", Config::Algo::DOUBLE_INDEXED}, {"1", Config::Algo::QUERY_INDEXED}, {"ctg", Config::Algo::CTG_SEED} };

const EMap<SequenceType> EnumTraits<SequenceType>::to_string = { {SequenceType::amino_acid,"prot"}, { SequenceType::nucleotide,"nucl"} };
const SEMap<SequenceType> EnumTraits<SequenceType>::from_string = { {"prot",SequenceType::amino_acid}, {"nucl",SequenceType::nucleotide} };

Config config;

pair<double, int> block_size(int64_t memory_limit, Sensitivity s, bool lin) {
	const double m = (double)memory_limit / 1e9;
	const int min = std::max(Search::sensitivity_traits[(int)align_mode.sequence_type].at(s).minimizer_window, 1),
		c = m < 40.0 && s <= Sensitivity::MORE_SENSITIVE && min == 1 ? 4 : 1;
	const double min_factor = std::min(1 / (double)min * 2, 1.0);
	const double max = s <= Sensitivity::DEFAULT ? 12.0 :
		(s <= Sensitivity::MORE_SENSITIVE ? 4.0 : 0.4);
	double b = m / (18.0 * min_factor / c + 2.0);
	/*if (b > 4)
		b = floor(b);
	else if (b > 0.4)
		b = floor(b * 10) / 10;
	else
		b = floor(b * 1000) / 1000;*/
	if (!config.no_block_size_limit && !lin)
		b = std::min(b, max);
	if (s >= Sensitivity::VERY_SENSITIVE)
		b = std::min(b, 2.1);
	//if (s >= Sensitivity::ULTRA_SENSITIVE)
		//b = std::min(b, 0.6);
	return { std::max(b, 0.001), c };
}

template<typename T>
T set_string_option(const string& s, const string& name, const vector<pair<string, T>>& values) {
	if (s.empty())
		return (T)0;
	for (auto& i : values)
		if (s == i.first)
			return i.second;
	throw std::runtime_error("Invalid argument for option " + name + ". Allowed values are:" + std::accumulate(values.begin(), values.end(), string(), [](const string& s, const pair<string, T>& v) { return s + ' ' + v.first; }));
}

void Config::set_sens(Sensitivity sens) {
	if (sensitivity != Sensitivity::DEFAULT)
		throw std::runtime_error("Sensitivity switches are mutually exclusive.");
	sensitivity = sens;
}

std::string Config::single_query_file() const {
	return query_file.empty() ? string() : query_file.front();
}

Compressor Config::compressor() const
{
	if (compression.empty() || compression == "0")
		return Compressor::NONE;
	else if (compression == "1")
		return Compressor::ZLIB;
	else if (compression == "zstd")
		return Compressor::ZSTD;
	else
		throw std::runtime_error("Invalid compression algorithm: " + compression);
}

Config::Config(int argc, const char **argv, bool check_io, CommandLineParser& parser)
{
   parser.add_command("makedb", "Build DIAMOND database from a FASTA file", makedb)
		.add_command("prepdb", "Prepare BLAST database for use with Diamond", prep_db)
		.add_command("blastp", "Align amino acid query sequences against a protein reference database", blastp)
		.add_command("blastx", "Align DNA query sequences against a protein reference database", blastx)
		.add_command("cluster", "Cluster protein sequences", cluster)
		.add_command("linclust", "Cluster protein sequences in linear time", LINCLUST)
		.add_command("realign", "Realign clustered sequences against their centroids", CLUSTER_REALIGN)
		.add_command("recluster", "Recompute clustering to fix errors", RECLUSTER)
		.add_command("reassign", "Reassign clustered sequences to the closest centroid", CLUSTER_REASSIGN)
		.add_command("view", "View DIAMOND alignment archive (DAA) formatted file", view)
		.add_command("merge-daa", "Merge DAA files", MERGE_DAA)
		.add_command("help", "Produce help message", help)
		.add_command("version", "Display version information", version)
		.add_command("getseq", "Retrieve sequences from a DIAMOND database file", getseq)
		.add_command("dbinfo", "Print information about a DIAMOND database file", dbinfo)
		.add_command("test", "Run regression tests", regression_test)
		.add_command("makeidx", "Make database index", makeidx)
		.add_command("greedy-vertex-cover", "Compute greedy vertex cover", GREEDY_VERTEX_COVER)
		.add_command("roc", "", roc)
		.add_command("benchmark", "", benchmark)
		.add_command("deepclust", "", DEEPCLUST)
#ifdef EXTRA
		.add_command("random-seqs", "", random_seqs)
		.add_command("sort", "", sort)
		.add_command("dbstat", "", db_stat)
		.add_command("mask", "", mask)
		.add_command("fastq2fasta", "", fastq2fasta)
		.add_command("read-sim", "", read_sim)
		.add_command("info", "", info)
		.add_command("seed-stat", "", seed_stat)
		.add_command("smith-waterman", "", smith_waterman)
		.add_command("translate", "", translate)
		.add_command("filter-blasttab", "", filter_blasttab)
		.add_command("show-cbs", "", show_cbs)
		.add_command("simulate-seqs", "", simulate_seqs)
		.add_command("split", "", split)
		.add_command("upgma", "", upgma)
		.add_command("upgmamc", "", upgma_mc)
		.add_command("reverse", "", reverse_seqs)
		.add_command("compute-medoids", "", compute_medoids)
		.add_command("mutate", "", mutate)
		.add_command("roc-id", "", rocid)
		.add_command("find-shapes", "", find_shapes)
		.add_command("composition", "", composition)
		.add_command("join", "", JOIN)
		.add_command("hashseqs", "", HASH_SEQS)
		.add_command("listseeds", "", LIST_SEEDS)
		.add_command("index-fasta", "", INDEX_FASTA)
		.add_command("fetch-seq", "", FETCH_SEQ)
		.add_command("blastn", "Align DNA query sequences against a DNA reference database", blastn)
		.add_command("length-sort", "", LENGTH_SORT)
		.add_command("wc", "", WORD_COUNT)
		.add_command("cut", "", CUT)
		.add_command("model-seqs", "", MODEL_SEQS)
#endif
		;

	auto& general = parser.add_group("General options", { makedb, blastp, blastx, cluster, view, prep_db, getseq, dbinfo, makeidx, CLUSTER_REALIGN, GREEDY_VERTEX_COVER, DEEPCLUST, RECLUSTER, MERGE_DAA, LINCLUST, CLUSTER_REASSIGN });
	general.add()
		("threads", 'p', "number of CPU threads", threads_)
		("verbose", 'v', "verbose console output", verbose)
		("log", 0, "enable debug log", debug_log)
		("quiet", 0, "disable console output", quiet)
		("tmpdir", 't', "directory for temporary files", tmpdir);

	auto& general_db = parser.add_group("General/database options", { makedb, blastp, blastx, cluster, prep_db, getseq, dbinfo, makeidx, CLUSTER_REALIGN, GREEDY_VERTEX_COVER, DEEPCLUST, RECLUSTER, LINCLUST, CLUSTER_REASSIGN });
	general_db.add()
		("db", 'd', "database file", database);

	auto& general_out = parser.add_group("General/output", { blastp, blastx, cluster, view, getseq, CLUSTER_REALIGN, GREEDY_VERTEX_COVER, DEEPCLUST, RECLUSTER, MERGE_DAA, LINCLUST, CLUSTER_REASSIGN });
	general_out.add()
		("out", 'o', "output file", output_file);

	auto& general_out2 = parser.add_group("General/output2", { blastp, blastx, cluster, view, CLUSTER_REALIGN, GREEDY_VERTEX_COVER, DEEPCLUST, RECLUSTER, LINCLUST, CLUSTER_REASSIGN });
	general_out2.add()
		("header", 0, "Use header lines in tabular output format (0/simple/verbose).", output_header, Option<vector<string>>(), 0);
	
    string dbstring;
	auto& makedb_opt = parser.add_group("Makedb options", { makedb, MERGE_DAA });
	makedb_opt.add()
		("in", 0, "input reference file in FASTA format/input DAA files for merge-daa", input_ref_file);

	auto& makedb_tax_opt = parser.add_group("Makedb/taxon options", { makedb });
	makedb_tax_opt.add()
		("taxonmap", 0, "protein accession to taxid mapping file", prot_accession2taxid)
		("taxonnodes", 0, "taxonomy nodes.dmp from NCBI", nodesdmp)
		("taxonnames", 0, "taxonomy names.dmp from NCBI", namesdmp);

	auto& align_clust_realign = parser.add_group("Aligner/Clustering/Realign options", { blastp, blastx, cluster, RECLUSTER, CLUSTER_REASSIGN, DEEPCLUST, CLUSTER_REALIGN, LINCLUST });
	align_clust_realign.add()
		("comp-based-stats", 0, "composition based statistics mode (0-4)", comp_based_stats, 1u)
		("masking", 0, "masking algorithm (none, seg, tantan=default)", masking_)
		("soft-masking", 0, "soft masking (none=default, seg, tantan)", soft_masking)
		("mmseqs-compat", 0, "", mmseqs_compat)
		("no-block-size-limit", 0, "", no_block_size_limit);

	auto& align_clust = parser.add_group("Aligner/Clustering options", { blastp, blastx, cluster, RECLUSTER, CLUSTER_REASSIGN, DEEPCLUST, LINCLUST });
	align_clust.add()		
		("evalue", 'e', "maximum e-value to report alignments (default=0.001)", max_evalue, 0.001)
		("motif-masking", 0, "softmask abundant motifs (0/1)", motif_masking)
		("approx-id", 0, "minimum approx. identity% to report an alignment/to cluster sequences", approx_min_id)
		("ext", 0, "Extension mode (banded-fast/banded-slow/full)", ext_);

	auto& aligner_view = parser.add_group("Aligner/view options", { blastp, blastx, view });
	aligner_view.add()
		("max-target-seqs", 'k', "maximum number of target sequences to report alignments for (default=25)", max_target_seqs_)
		("top", 0, "report alignments within this percentage range of top alignment score (overrides --max-target-seqs)", toppercent, 100.0);

	auto& aligner_sens = parser.add_group("Aligner/sens options", { blastp, blastx, makeidx });
	aligner_sens.add()
		("faster", 0, "enable faster mode", mode_faster)
		("fast", 0, "enable fast mode", mode_fast)
		("mid-sensitive", 0, "enable mid-sensitive mode", mode_mid_sensitive)
		("sensitive", 0, "enable sensitive mode)", mode_sensitive)
		("more-sensitive", 0, "enable more sensitive mode", mode_more_sensitive)
		("very-sensitive", 0, "enable very sensitive mode", mode_very_sensitive)
		("ultra-sensitive", 0, "enable ultra sensitive mode", mode_ultra_sensitive)
		("shapes", 's', "number of seed shapes (default=all available)", shapes);

	auto& aligner = parser.add_group("Aligner options", { blastp, blastx });
	aligner.add()
		("query", 'q', "input query file", query_file)
		("strand", 0, "query strands to search (both/minus/plus)", query_strands, string("both"))
		("un", 0, "file for unaligned queries", unaligned)
		("al", 0, "file or aligned queries", aligned_file)
		("unfmt", 0, "format of unaligned query file (fasta/fastq)", unfmt, string("fasta"))
		("alfmt", 0, "format of aligned query file (fasta/fastq)", alfmt, string("fasta"))
		("unal", 0, "report unaligned queries (0=no, 1=yes)", report_unaligned, -1)
		("max-hsps", 0, "maximum number of HSPs per target sequence to report for each query (default=1)", max_hsps, 1u)
		("range-culling", 0, "restrict hit culling to overlapping query ranges", query_range_culling)
		("compress", 0, "compression for output files (0=none, 1=gzip, zstd)", compression)
		("min-score", 0, "minimum bit score to report alignments (overrides e-value setting)", min_bit_score)
		("id", 0, "minimum identity% to report an alignment", min_id)
		("query-cover", 0, "minimum query cover% to report an alignment", query_cover)
		("subject-cover", 0, "minimum subject cover% to report an alignment", subject_cover)
		("swipe", 0, "exhaustive alignment against all database sequences", swipe_all)
		("iterate", 0, "iterated search with increasing sensitivity", iterate, Option<vector<string>>(), 0)
		("global-ranking", 'g', "number of targets for global ranking", global_ranking_targets)
		("block-size", 'b', "sequence block size in billions of letters (default=2.0)", chunk_size)
		("index-chunks", 'c', "number of chunks for index processing (default=4)", lowmem_)
		("parallel-tmpdir", 0, "directory for temporary files used by multiprocessing", parallel_tmpdir)
		("gapopen", 0, "gap open penalty", gap_open, -1)
		("gapextend", 0, "gap extension penalty", gap_extend, -1)
		("matrix", 0, "score matrix for protein alignment (default=BLOSUM62)", matrix, string("blosum62"))
		("custom-matrix", 0, "file containing custom scoring matrix", matrix_file)
		("frameshift", 'F', "frame shift penalty (default=disabled)", frame_shift)
		("long-reads", 0, "short for --range-culling --top 10 -F 15", long_reads)
		("query-gencode", 0, "genetic code to use to translate query (see user manual)", query_gencode, 1u)
		("salltitles", 0, "include full subject titles in DAA file", salltitles)
		("sallseqid", 0, "include all subject ids in DAA file", sallseqid)
		("no-self-hits", 0, "suppress reporting of identical self hits", no_self_hits)
		("taxonlist", 0, "restrict search to list of taxon ids (comma-separated)", taxonlist)
		("taxon-exclude", 0, "exclude list of taxon ids (comma-separated)", taxon_exclude)
		("seqidlist", 0, "filter the database by list of accessions", seqidlist)
		("skip-missing-seqids", 0, "ignore accessions missing in the database", skip_missing_seqids);

	auto& format = parser.add_group("Output format options", { blastp, blastx, view, CLUSTER_REALIGN });
	format.add()
		("outfmt", 'f', "output format\n\
\t0   = BLAST pairwise\n\
\t5   = BLAST XML\n\
\t6   = BLAST tabular\n\
\t100 = DIAMOND alignment archive (DAA)\n\
\t101 = SAM\n\
\t102 = Taxonomic classification\n\
\t103 = PAF\n\
\t104 = JSON (flat)\n\n\
\tValues 6 and 104 may be followed by a space-separated list of these keywords:\n\n\
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
\tqseq_translated means Aligned part of query sequence (translated)\n\
\tfull_qseq means Query sequence\n\
\tfull_qseq_mate means Query sequence of the mate\n\
\tsseq means Aligned part of subject sequence\n\
\tfull_sseq means Subject sequence\n\
\tevalue means Expect value\n\
\tbitscore means Bit score\n\
\tcorrected_bitscore means Bit score corrected for edge effects\n\
\tscore means Raw score\n\
\tlength means Alignment length\n\
\tpident means Percentage of identical matches\n\
\tapprox_pident means Approximate percentage of identical matches\n\
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
\n\tDefault: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", output_format);

	auto& cluster_opt = parser.add_group("Clustering options", { cluster, RECLUSTER, DEEPCLUST, LINCLUST });
	kmer_ranking = false;
	cluster_opt.add()
		("cluster-steps", 0, "Clustering steps", cluster_steps)
#ifdef KEEP_TARGET_ID
		("kmer-ranking", 0, "Rank sequences based on kmer frequency in linear stage", kmer_ranking)
#endif
		;

	auto& cluster_reassign_opt = parser.add_group("Clustering/reassign options", { cluster, RECLUSTER, CLUSTER_REASSIGN, GREEDY_VERTEX_COVER, DEEPCLUST, LINCLUST });
	cluster_reassign_opt.add()
		("memory-limit", 'M', "Memory limit in GB (default = 16G)", memory_limit)
		("member-cover", 0, "Minimum coverage% of the cluster member sequence (default=80.0)", member_cover, 80.0);

	auto& gvc_opt = parser.add_group("GVC options", { GREEDY_VERTEX_COVER });
	gvc_opt.add()
		("centroid-out", 0, "Output file for centroids", centroid_out)
		("edges", 0, "Input file for greedy vertex cover", edges)
		("edge-format", 0, "Edge format for greedy vertex cover (default/triplet)", edge_format);

	auto& realign_opt = parser.add_group("Cluster input options", { CLUSTER_REALIGN, RECLUSTER, CLUSTER_REASSIGN });
	realign_opt.add()
		("clusters", 0, "Clustering input file mapping sequences to representatives", clustering);

	string algo_str;

	auto& advanced_gen = parser.add_group("Advanced/general", { blastp, blastx, blastn, CLUSTER_REASSIGN, regression_test, cluster, DEEPCLUST, LINCLUST, makedb });
	advanced_gen.add()
		("file-buffer-size", 0, "file buffer size in bytes (default=67108864)", file_buffer_size, (size_t)67108864)
		("no-unlink", 0, "Do not unlink temporary files.", no_unlink)
		("ignore-warnings", 0, "Ignore warnings", ignore_warnings)
		("no-parse-seqids", 0, "Print raw seqids without parsing", no_parse_seqids);

	auto& advanced_aln_cluster = parser.add_group("Advanced options aln/cluster", { blastp, blastx, blastn, CLUSTER_REASSIGN, regression_test, cluster, DEEPCLUST, LINCLUST, RECLUSTER });
	advanced_aln_cluster.add()
		("bin", 0, "number of query bins for seed search", query_bins_)
		("ext-chunk-size", 0, "chunk size for adaptive ranking (default=auto)", ext_chunk_size)
		("no-ranking", 0, "disable ranking heuristic", no_ranking)
		("dbsize", 0, "effective database size (in letters)", db_size)
		("no-auto-append", 0, "disable auto appending of DAA and DMND file extensions", no_auto_append)
		("tantan-minMaskProb", 0, "minimum repeat probability for masking (default=0.9)", tantan_minMaskProb, 0.9);

	auto& advanced = parser.add_group("Advanced options", { blastp, blastx, blastn, regression_test });
	advanced.add()
		("algo", 0, "Seed search algorithm (0=double-indexed/1=query-indexed/ctg=contiguous-seed)", algo_str)
		("min-orf", 'l', "ignore translated sequences without an open reading frame of at least this length", run_len)
		("seed-cut", 0, "cutoff for seed complexity", seed_cut_)
		("freq-masking", 0, "mask seeds based on frequency", freq_masking)
		("freq-sd", 0, "number of standard deviations for ignoring frequent seeds", freq_sd_, 0.0)
		("id2", 0, "minimum number of identities for stage 1 hit", min_identities_)
		("linsearch", 0, "only consider seed hits against longest target for identical seeds", linsearch)
		("lin-stage1", 0, "only consider seed hits against longest query for identical seeds", lin_stage1)
		("xdrop", 'x', "xdrop for ungapped alignment", ungapped_xdrop, 12.3)
		("gapped-filter-evalue", 0, "E-value threshold for gapped filter (auto)", gapped_filter_evalue_, -1.0)
		("band", 0, "band for dynamic programming computation", padding)
		("shape-mask", 0, "seed shapes", shape_mask)
		("multiprocessing", 0, "enable distributed-memory parallel processing", multiprocessing)
		("mp-init", 0, "initialize multiprocessing run", mp_init)
		("mp-recover", 0, "enable continuation of interrupted multiprocessing run", mp_recover)
		("mp-query-chunk", 0, "process only a single query chunk as specified", mp_query_chunk, -1)
		("culling-overlap", 0, "minimum range overlap with higher scoring hit to delete a hit (default=50%)", inner_culling_overlap, 50.0)
		("taxon-k", 0, "maximum number of targets to report per species", taxon_k, (uint64_t)0)
		("range-cover", 0, "percentage of query range to be covered for range culling (default=50%)", query_range_cover, 50.0)
		("xml-blord-format", 0, "Use gnl|BL_ORD_ID| style format in XML output", xml_blord_format)
		("sam-query-len", 0, "add the query length to the SAM format (tag ZQ)", sam_qlen_field)
		("stop-match-score", 0, "Set the match score of stop codons against each other.", stop_match_score, 1)		
		("target-indexed", 0, "Enable target-indexed mode", target_indexed)
		("unaligned-targets", 0, "", unaligned_targets)
		("cut-bar", 0, "", cut_bar)
		("check-multi-target", 0, "", check_multi_target)
		("roc-file", 0, "", roc_file)
		("family-map", 0, "", family_map)
		("family-map-query", 0, "", family_map_query)
		("query-parallel-limit", 0, "", query_parallel_limit, 3000000u)
		("log-evalue-scale", 0, "", log_evalue_scale, 1.0 / std::log(2.0))
		("bootstrap", 0, "", bootstrap)
		("heartbeat", 0, "", heartbeat)
		("mp-self", 0, "", mp_self)
#ifdef EXTRA
            ("zdrop", 'z', "zdrop for gapped dna alignment", zdrop, 40)
#endif
		("query-or-subject-cover", 0, "", query_or_target_cover);

	auto& view_align_options = parser.add_group("View/Align options", { view, blastp, blastx });
	view_align_options.add()
		("daa", 'a', "DIAMOND alignment archive (DAA) file", daa_file);

	auto& view_options = parser.add_group("View options", { view });
	view_options.add()
		("forwardonly", 0, "only show alignments of forward strand", forwardonly);

	auto& getseq_options = parser.add_group("Getseq options", { getseq });
	getseq_options.add()
		("seq", 0, "Space-separated list of sequence numbers to display.", seq_no);

	double rank_ratio2, lambda, K;
	unsigned window, min_ungapped_score, hit_band, min_hit_score;
	auto& deprecated_options = parser.add_group("", { blastp, blastx });
	deprecated_options.add()
		("window", 'w', "window size for local hit search", window)
		("ungapped-score", 0, "minimum alignment score to continue local extension", min_ungapped_score)
		("hit-band", 0, "band for hit verification", hit_band)
		("hit-score", 0, "minimum score to keep a tentative alignment", min_hit_score)
		("gapped-xdrop", 'X', "xdrop for gapped alignment in bits", gapped_xdrop, 20.0)
		("rank-ratio2", 0, "include subjects within this ratio of last hit (stage 2)", rank_ratio2, -1.0)
		("rank-ratio", 0, "include subjects within this ratio of last hit", rank_ratio, -1.0)
		("lambda", 0, "lambda parameter for custom matrix", lambda)
		("K", 0, "K parameter for custom matrix", K);

	double query_match_distance_threshold;
	double length_ratio_threshold;
	double cbs_angle;

#ifdef EXTRA
	auto& hidden_options = parser.add_group("", {});
#else
	auto& hidden_options = parser.add_group("", {}, true);
#endif
	hidden_options.add()
		("match1", 0, "", match_file1)
		("match2", 0, "", match_file2)
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
		("path-cutoff", 0, "", path_cutoff, 0.92)
		("sw", 0, "", use_smith_waterman)
		("superblock", 0, "", superblock, 128)
		("max-cells", 0, "", max_cells, 10000000u)
		("load-balancing", 0, "", load_balancing, (unsigned)Config::query_parallel)
		("log-query", 0, "", log_query)
		("log-subject", 0, "", log_subject)
		("palign", 0, "", threads_align)
		("score-ratio", 0, "", score_ratio, 0.9)
		("fetch-size", 0, "trace point fetch size", fetch_size, 4096u)
		("target-fetch-size", 0, "number of target sequences to fetch for seed extension", target_fetch_size, 4u)
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
		("use-dataset-field", 0, "", use_dataset_field)
		("store-query-quality", 0, "", store_query_quality)
		("swipe-chunk-size", 0, "", swipe_chunk_size, 256u)
		("hard-masked", 0, "", hardmasked)
		("cbs-window", 0, "", cbs_window, 40)
		("no-dict", 0, "", no_dict)
		("upgma-edge-limit", 0, "", upgma_edge_limit, (uint64_t)10000000)
		("tree", 0, "", tree_file)
		("upgma-dist", 0, "", upgma_dist)
		("upgma-input", 0, "", upgma_input)
		("log-extend", 0, "", log_extend)
		("chaining-maxgap", 0, "", chaining_maxgap, 2000)
		("tantan-maxRepeatOffset", 0, "maximum tandem repeat period to consider (50)", tantan_maxRepeatOffset, 15)
		("tantan-ungapped", 0, "use tantan masking in ungapped mode", tantan_ungapped)
		("chaining-range-cover", 0, "", chaining_range_cover, (size_t)8)
		("no-swipe-realign", 0, "", no_swipe_realign)
		("chaining-maxnodes", 0, "", chaining_maxnodes)
		("cutoff-score-8bit", 0, "", cutoff_score_8bit, 240)
		("min-band-overlap", 0, "", min_band_overlap, 0.0)
		("min-realign-overhang", 0, "", min_realign_overhang, 30)
		("ungapped-window", 0, "", ungapped_window, 48)
		("gapped-filter-diag-score", 0, "", gapped_filter_diag_bit_score, 12.0)
		("gapped-filter-window", 0, "", gapped_filter_window, 200)
		("output-hits", 0, "", output_hits)
		("ungapped-evalue", 0, "", ungapped_evalue_, -1.0)
		("ungapped-evalue-short", 0, "", ungapped_evalue_short_, -1.0)
		("no-logfile", 0, "", no_logfile)
		("band-bin", 0, "", band_bin, 24)
		("col-bin", 0, "", col_bin, 400)
		("self", 0, "", self)
		("trace-pt-fetch-size", 0, "", trace_pt_fetch_size, (int64_t)10e9)
		("tile-size", 0, "", tile_size, (uint32_t)1024)
		("short-query-ungapped-bitscore", 0, "", short_query_ungapped_bitscore, 25.0)
		("short-query-max-len", 0, "", short_query_max_len, 60)
		("gapped-filter-evalue1", 0, "", gapped_filter_evalue1, 2000.0)
		("ext-yield", 0, "", ext_min_yield)
		("full-sw-len", 0, "", full_sw_len)
		("relaxed-evalue-factor", 0, "", relaxed_evalue_factor, 1.0)
		("type", 0, "", type)
		("raw", 0, "", raw)
		("chaining-len-cap", 0, "", chaining_len_cap, 2.0)
		("chaining-min-nodes", 0, "", chaining_min_nodes, (size_t)200)
		("fast-tsv", 0, "", fast_tsv)
		("target-parallel-verbosity", 0, "", target_parallel_verbosity, UINT_MAX)
		("query-memory", 0, "", query_memory)
		("memory-intervals", 0, "", memory_intervals, (size_t)2)
		("seed-hit-density", 0, "", seedhit_density)
		("chunk-size-multiplier", 0, "", chunk_size_multiplier, (size_t)4)
		("score-drop-factor", 0, "", ranking_score_drop_factor, 0.95)
		("left-most-interval", 0, "", left_most_interval, 32)
		("ranking-cutoff-bitscore", 0, "", ranking_cutoff_bitscore, 25.0)
		("no-forward-fp", 0, "", no_forward_fp)
		("no-ref-masking", 0, "", no_ref_masking)
		("target-bias", 0, "", target_bias)
		("output-fp", 0, "", output_fp)
		("family-cap", 0, "", family_cap)
		("cbs-matrix-scale", 0, "", cbs_matrix_scale, 1)
		("query-count", 0, "", query_count, (size_t)1)
		("cbs-angle", 0, "", cbs_angle, -1.0)
		("cbs-err-tolerance", 0, "", cbs_err_tolerance, 0.00000001)
		("cbs-it-limit", 0, "", cbs_it_limit, 2000)
		("hash_join_swap", 0, "", hash_join_swap)
		("deque_bucket_size", 0, "", deque_bucket_size, (size_t)524288)
		("query-match-distance-threshold", 0, "", query_match_distance_threshold, -1.0)
		("length-ratio-threshold", 0, "", length_ratio_threshold, -1.0)
		("max-swipe-dp", 0, "", max_swipe_dp, (int64_t)1000000)
		("short-seqids", 0, "", short_seqids)
		("no-reextend", 0, "", no_reextend)
		("no-reorder", 0, "", no_reorder)
		("file1", 0, "", file1)
		("file2", 0, "", file2)
		("key2", 0, "", key2)
		("motif-mask-file", 0, "", motif_mask_file)
		("max-motif-len", 0, "", max_motif_len, 30)
		("chaining-stacked-hsp-ratio", 0, "", chaining_stacked_hsp_ratio, 0.5)
		("swipe-task-size", 0, "", swipe_task_size, (int64_t)100000000)
		("minimizer-window", 0, "", minimizer_window_)
		("min_task_trace_pts", 0, "", min_task_trace_pts, (int64_t)1024)
		("sketch-size", 0, "", sketch_size)
		("oid-list", 0, "", oid_list)
		("bootstrap-block", 0, "", bootstrap_block, (int64_t)1000000)
		("centroid-factor", 0, "", centroid_factor, (int64_t)3)
		("timeout", 0, "", timeout)
		("resume", 0, "", resume)
		("target_hard_cap", 0, "", target_hard_cap)
		("mapany", 0, "", mapany)
		("neighbors", 0, "", neighbors)
		("reassign-overlap", 0, "", reassign_overlap, 0.3)
		("reassign-ratio", 0, "", reassign_ratio, 0.5)
		("reassign-max", 0, "", reassign_max)
		("add-self-aln", 0, "", add_self_aln)
		("weighted-gvc", 0, "", weighted_gvc)
		("hamming-ext", 0, "", hamming_ext)
		("diag-filter-id", 0, "", diag_filter_id)
		("diag-filter-cov", 0, "", diag_filter_cov)
		("strict-gvc", 0, "", strict_gvc)
		("dbtype", 0, "type of sequences in database file (nucl/prot)", dbstring, string("prot"))
		("penalty", 0, "blastn mismatch penalty", mismatch_penalty, -3)
		("reward", 0, "blastn match reward", match_reward, 2)
		("cluster-similarity", 0, "Clustering similarity measure (default=\"normalized_bitscore_global\")", cluster_similarity)
		("cluster-threshold", 0, "Threshold for the similarity measure (default=50%)", cluster_threshold)
		("cluster-graph-file", 0, "Filename for dumping the graph or reading the graph if cluster-restart", cluster_graph_file)
		("cluster-restart", 0, "Restart clustering from dumped graph", cluster_restart)
		("mcl-expansion", 0, "MCL expansion coefficient (default=2)", cluster_mcl_expansion, (uint32_t)2)
		("mcl-inflation", 0, "MCL inflation coefficient (default=2.0)", cluster_mcl_inflation, 2.0)
		("mcl-chunk-size", 0, "MCL chunk size per thread (default=100)", cluster_mcl_chunk_size, (uint32_t)1)
		("mcl-max-iterations", 0, "MCL maximum iterations (default=100)", cluster_mcl_max_iter, (uint32_t)100)
		("mcl-sparsity-switch", 0, "MCL switch to sparse matrix computation (default=0.8) ", cluster_mcl_sparsity_switch, 0.8)
		("mcl-nonsymmetric", 0, "Do not symmetrize the transistion matrix before clustering", cluster_mcl_nonsymmetric)
		("mcl-stats", 0, "Some stats about the connected components in MCL", cluster_mcl_stats)
		("cluster-algo", 0, "Clustering algorithm (\"mcl\")", cluster_algo)
		("approx-backtrace", 0, "", approx_backtrace)
		("prefix-scan", 0, "", prefix_scan)
		("narrow-band-cov", 0, "", narrow_band_cov)
		("narrow-band-factor", 0, "", narrow_band_factor)
		("anchor-window", 0, "", anchor_window, 12)
		("anchor-score", 0, "", anchor_score, 1.0)
		("classic-band", 0, "", classic_band)
		("no_8bit_extension", 0, "", no_8bit_extension)
		("anchored-swipe", 0, "", anchored_swipe)
		("no_chaining_merge_hsps", 0, "", no_chaining_merge_hsps)
		("recluster_bd", 0, "", recluster_bd)
		("pipeline-short", 0, "", pipeline_short)
		("graph-algo", 0, "", graph_algo, string("gvc"))
		("tsv-read-size", 0, "", tsv_read_size, int64_t(GIGABYTES))
#ifndef KEEP_TARGET_ID
		("kmer-ranking", 0, "Rank sequences based on kmer frequency in linear stage", kmer_ranking)
#endif
		;

	parser.store(argc, argv, command);

	if (debug_log)
		verbosity = 3;
	else if (quiet)
		verbosity = 0;
	else if (verbose)
		verbosity = 2;
	else if (((command == Config::view || command == blastx || command == blastp || command == blastn) && output_file == "" && argc != 2)
		|| command == Config::version || command == getseq || command == fastq2fasta || command == regression_test)
		verbosity = 0;
	else
		verbosity = 1;

	if (verbosity >= 1 || command == regression_test) {
		ostream& header_out = command == Config::help ? cout : cerr;
		header_out << Const::program_name << " v" << Const::version_string << "." << (unsigned)Const::build_version << " (C) Max Planck Society for the Advancement of Science, Benjamin Buchfink, University of Tuebingen" << endl;
		header_out << "Documentation, support and updates available at http://www.diamondsearch.org" << endl;
		header_out << "Please cite: http://dx.doi.org/10.1038/s41592-021-01101-x Nature Methods (2021)" << endl << endl;
	}
	log_stream << Const::program_name << " v" << Const::version_string << "." << (unsigned)Const::build_version << endl;

    if(argc == 2 && command != version && command != regression_test) {
        if (command != help) {
            parser.print_documentation(command);
        } else {
            parser.print_help();
        }
    }


	if (toppercent != 100.0 && max_target_seqs_.present())
		throw std::runtime_error("--top and --max-target-seqs are mutually exclusive.");

	if (command == blastx && no_self_hits)
		throw std::runtime_error("--no-self-hits option is not supported in blastx mode.");

	if (long_reads) {
		query_range_culling = true;
		if (toppercent == 100.0)
			toppercent = 10.0;
		if (frame_shift == 0)
			frame_shift = 15;
	}

	if (global_ranking_targets > 0 && (query_range_culling || taxon_k || multiprocessing || mp_init || mp_recover || comp_based_stats >= 2 || frame_shift > 0))
		throw std::runtime_error("Global ranking is not supported in this mode.");



#ifdef EXTRA
	if (comp_based_stats >= Stats::CBS::COUNT)
#else
	if (comp_based_stats >= 5)
#endif
		throw std::runtime_error("Invalid value for --comp-based-stats. Permitted values: 0, 1, 2, 3, 4.");

	Stats::comp_based_stats = Stats::CBS(comp_based_stats, query_match_distance_threshold, length_ratio_threshold, cbs_angle);

	if (command == blastx && !Stats::CBS::support_translated(comp_based_stats))
		throw std::runtime_error("This mode of composition based stats is not supported for translated searches.");

    if (check_io) {
		switch (command) {
		case Config::makedb: {
            if (database == "")
                throw std::runtime_error("Missing parameter: database file (--db/-d)");
            if (chunk_size != 0.0)
                throw std::runtime_error(
                        "Invalid option: --block-size/-b. Block size is set for the alignment commands.");
            break;
        }
		case Config::blastp:
		case Config::blastx:
        case Config::blastn:
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
				if (!compression.empty())
					throw std::runtime_error("Compression is not supported for DAA format.");
				if (!no_auto_append)
					auto_append_extension(output_file, ".daa");
			}
			break;
		default:
			;
		}

		switch (command) {
		case Config::dbinfo:
			if (database == "")
				throw std::runtime_error("Missing parameter: database file (--db/-d)");
		}
	}

	switch (verbosity) {
	case 0:
		message_stream = MessageStream(false);
		break;
	case 3:
		log_stream = MessageStream(true, !config.no_logfile);
		verbose_stream = MessageStream(true, !config.no_logfile);
		message_stream = MessageStream(true, !config.no_logfile);
		break;
	case 2:
		verbose_stream = MessageStream();
	default:
		;
	}

	invocation = join(" ", vector<string>(&argv[0], &argv[argc]));
	log_stream << invocation << endl;

	if (!no_auto_append) {
		if (command == Config::makedb)
			auto_append_extension(database, ".dmnd");
		if (command == Config::view)
			auto_append_extension(daa_file, ".daa");
		if (compression == "1")
			auto_append_extension(output_file, ".gz");
		if (compression == "zstd")
			auto_append_extension(output_file, ".zst");
	}

	
#ifndef NDEBUG
	verbose_stream << "Assertions enabled." << endl;
#endif
	set_option(threads_, (int)std::thread::hardware_concurrency());

	switch (command) {
	case Config::makedb:
	case Config::blastp:
	case Config::blastx:
    case Config::blastn:
	case Config::view:
	case Config::cluster:
	case Config::DEEPCLUST:
	case Config::LINCLUST:
	case Config::regression_test:
	case Config::compute_medoids:
	case Config::CLUSTER_REASSIGN:
	case Config::GREEDY_VERTEX_COVER:
	case Config::RECLUSTER:
		if (argc != 2)
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
	case Config::DEEPCLUST:
	case Config::LINCLUST:
	case Config::regression_test:
	case Config::compute_medoids:
	case Config::LIST_SEEDS:
	case Config::CLUSTER_REASSIGN:
	case Config::CLUSTER_REALIGN:
	case Config::RECLUSTER:
	case Config::MODEL_SEQS:
		if (frame_shift != 0 && command == Config::blastp)
			throw std::runtime_error("Frameshift alignments are only supported for translated searches.");
		if (query_range_culling && frame_shift == 0)
			throw std::runtime_error("Query range culling is only supported in frameshift alignment mode (option -F).");
		if (matrix_file == "") {
			score_matrix = ScoreMatrix(to_upper_case(matrix), gap_open, gap_extend, frame_shift, stop_match_score, 0, cbs_matrix_scale, mmseqs_compat);
		}
		else {
			if (gap_open == -1 || gap_extend == -1)
				throw std::runtime_error("Custom scoring matrices require setting the --gapopen and --gapextend options.");
			if (!output_format.empty() && (output_format.front() == "daa" || output_format.front() == "100"))
				throw std::runtime_error("Custom scoring matrices are not supported for the DAA format.");
			if (comp_based_stats > 1)
				throw std::runtime_error("This value for --comp-based-stats is not supported when using a custom scoring matrix.");
			score_matrix = ScoreMatrix(matrix_file, gap_open, gap_extend, stop_match_score, ScoreMatrix::Custom());
		}
		if (argc != 2)
			message_stream << "Scoring parameters: " << score_matrix << endl;
		Masking::instance = unique_ptr<Masking>(new Masking(score_matrix));
	}

	if (command == Config::blastp || command == Config::blastx || command == Config::blastn || command == Config::benchmark || command == Config::model_sim || command == Config::opt
		|| command == Config::mask || command == Config::cluster || command == Config::compute_medoids || command == Config::regression_test || command == Config::CLUSTER_REASSIGN
		|| command == Config::RECLUSTER || command == Config::DEEPCLUST || command == Config::LINCLUST) {
		if (tmpdir == "")
			tmpdir = extract_dir(output_file);

		raw_ungapped_xdrop = score_matrix.rawscore(ungapped_xdrop);
		verbose_stream << "CPU features detected: " << SIMD::features() << endl;
		log_stream << "L3 cache size: " << l3_cache_size() << endl;
	}

	sensitivity = Sensitivity::DEFAULT;
	if (mode_faster) set_sens(Sensitivity::FASTER);
	if (mode_fast) set_sens(Sensitivity::FAST);
	if (mode_mid_sensitive) set_sens(Sensitivity::MID_SENSITIVE);
	if (mode_sensitive) set_sens(Sensitivity::SENSITIVE);
	if (mode_more_sensitive) set_sens(Sensitivity::MORE_SENSITIVE);
	if (mode_very_sensitive) set_sens(Sensitivity::VERY_SENSITIVE);
	if (mode_ultra_sensitive) set_sens(Sensitivity::ULTRA_SENSITIVE);

	algo = from_string<Algo>(algo_str);
    dbtype = from_string<SequenceType>(dbstring);
	Translator::init(query_gencode);

	if (command == blastx || command == blastn)
		input_value_traits = nucleotide_traits;

	if (query_strands != "both" && query_strands != "minus" && query_strands != "plus")
		throw std::runtime_error("Invalid value for parameter --strand");

	if (unfmt == "fastq" || alfmt == "fastq")
		store_query_quality = true;
	if (!aligned_file.empty())
		log_stream << "Aligned file format: " << alfmt << endl;

	if (command == blastx) {
		if (query_file.size() > 2)
			throw std::runtime_error("A maximum of 2 query files is supported in blastx mode.");
	}
	else if (query_file.size() > 1)
		throw std::runtime_error("--query/-q has more than one argument.");

	if (target_indexed && lowmem_ != 1)
		throw std::runtime_error("--target-indexed requires -c1.");

	/*log_stream << "sizeof(hit)=" << sizeof(hit) << " sizeof(packed_uint40_t)=" << sizeof(packed_uint40_t)
		<< " sizeof(sorted_list::entry)=" << sizeof(sorted_list::entry) << endl;*/

	if (swipe_all) {
		algo = Algo::DOUBLE_INDEXED;
	}

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
#ifdef KEEP_TARGET_ID
	log_stream << " KEEP_TARGET_ID";
#endif
	log_stream << endl;
}