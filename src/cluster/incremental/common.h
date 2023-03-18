#pragma once
#include "../../data/sequence_file.h"
#include "../search/search.h"
#include "../../data/fasta/fasta_file.h"
#include "../util/tsv/file.h"

namespace Cluster { namespace Incremental {

struct Config {
	Config();
	~Config(); 
	MessageStream                       message_stream;
	int                                 verbosity;
	std::vector<Sensitivity>            sens;
	double                              block_size;
	std::unique_ptr<OutputFormat>       output_format;
	std::unique_ptr<SequenceFile>       db;
	std::shared_ptr<FastaFile>          centroids;
	std::unique_ptr<Util::Tsv::File>    output_file;
	TaskTimer                          total_time;
	int64_t                             seqs_processed;
	int64_t                             letters_processed;
	std::vector<CentroidId>             oid2centroid;
	std::vector<OId>                    centroid2oid;
	std::vector<std::unique_ptr<Block>> cache;
	int64_t                             time_self_aln;
	std::vector<int64_t>                time_search;
	std::vector<int64_t>                problem_size;
	int64_t                             problem_size_self;

	void status_msg();
	void save_state();
	void load_state();
};

}}