#include <iostream>
#include <vector>
#include "tsv_record.h"
#include "../basic/config.h"
#include "../basic/value.h"
#include "../util/util.h"
#include "../basic/sequence.h"
#include "../util/seq_file_format.h"
#include "../util/sequence/sequence.h"
#include "../stats/cbs.h"

using std::cout;
using std::endl;
using std::vector;

void filter_blasttab() {
	TextInputFile in("");
	TSVRecord r;
	string query;
	size_t query_hit;
	while (in >> r) {
		if (r.qseqid != query) {
			query = r.qseqid;
			query_hit = 0;
		}
		else
			++query_hit;
		if(query_hit < config.max_alignments && r.evalue <= config.max_evalue)
			cout << r << endl;
	}
}

void split() {
	TextInputFile in(config.single_query_file());
	string id;
	vector<Letter> seq;
	size_t n = 0, f = 0, b = (size_t)(config.chunk_size * 1e9);
	OutputFile *out = new OutputFile(std::to_string(f) + ".faa.gz", Compressor::ZLIB);
	while (FASTA_format().get_seq(id, seq, in, value_traits)) {
		if (n >= b) {
			out->close();
			delete out;
			out = new OutputFile(std::to_string(++f) + ".faa.gz", Compressor::ZLIB);
			n = 0;
		}
		string blast_id = Util::Seq::seqid(id.c_str(), false);
		Util::Seq::format(Sequence(seq), blast_id.c_str(), nullptr, *out, "fasta", amino_acid_traits);
		n += seq.size();
	}
	out->close();
	delete out;
}

void composition() {
	TextInputFile in(config.single_query_file());
	string id;
	vector<Letter> seq;
	while (FASTA_format().get_seq(id, seq, in, value_traits)) {
		auto c = Stats::composition(seq);
		for (double x : c)
			std::cout << x << '\t';
		std::cout << endl;
	}
}