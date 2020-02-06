#include <iostream>
#include <string>
#include <limits.h>
#include "../util/io/temp_file.h"
#include "../util/io/text_input_file.h"
#include "test.h"
#include "../util/sequence/sequence.h"
#include "../util/log_stream.h"
#include "../data/reference.h"
#include "../basic/statistics.h"
#include "../run/workflow.h"

using std::endl;
using std::string;

namespace Test {

void run() {
	task_timer timer("Generating test dataset");
	TempFile proteins;
	for (size_t i = 0; i < sizeof(seqs)/sizeof(seqs[0]); ++i)
		Util::Sequence::format(sequence::from_string(seqs[i][1]), seqs[i][0], nullptr, proteins, "fasta", amino_acid_traits);
	TextInputFile query_file(proteins);
	timer.finish();

	const char* args[] = { "diamond", "blastp" };
	config = Config(2, args, false);
	config.command = Config::makedb;
	TempFile *db_file;
	make_db(&db_file, &query_file);
	DatabaseFile db(*db_file);

	statistics.reset();
	config.command = Config::blastp;
	config.algo = 0;
	config.max_alignments = std::numeric_limits<uint64_t>::max();
	config.mode_more_sensitive = true;
	config.lowmem = 1;
	
	Workflow::Search::Options opt;
	opt.db = &db;
	query_file.rewind();
	opt.query_file = &query_file;
	TempFile output_file;
	opt.consumer = &output_file;

	Workflow::Search::run(opt);

	InputFile out_in(output_file);
	cout << out_in.hash() << endl;

	out_in.close_and_delete();
	query_file.close_and_delete();
	db.close();
	delete db_file;
}

}