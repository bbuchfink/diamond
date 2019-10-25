#include <vector>
#include "tsv_record.h"
#include "../basic/config.h"
#include "../util/math/sparse_matrix.h"
#include "../basic/value.h"
#include "../util/util.h"
#include "../basic/sequence.h"
#include "../util/seq_file_format.h"
#include "../util/sequence/sequence.h"

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

void mcl() {
	SparseMatrix graph(std::cin);
	graph.print_stats();
}

vector<char> generate_random_seq(size_t length)
{
	vector<char> seq;
	seq.reserve(length);
	for (size_t i = 0; i < length; ++i)
		seq.push_back(get_distribution<20>(background_freq));
	return seq;
}

void simulate_seqs() {
	const size_t l = 300, n = atoi(config.seq_no[0].c_str());
	for (size_t i = 0; i < n; ++i) {
		cout << ">" << i << endl;
		cout << sequence(generate_random_seq(l)) << endl;
	}
}

void split() {
	TextInputFile in(config.query_file);
	vector<char> id, seq;
	size_t n = 0, f = 0, b = (size_t)(config.chunk_size * 1e9);
	OutputFile *out = new OutputFile(std::to_string(f) + ".faa.gz", true);
	while (FASTA_format().get_seq(id, seq, in)) {
		if (seq.size() + n > b) {
			out->close();
			delete out;
			out = new OutputFile(std::to_string(++f) + ".faa.gz", true);
			n = 0;
		}
		string blast_id = ::blast_id(string(id.data(), id.size()));
		Util::Sequence::format(sequence(seq), blast_id.c_str(), nullptr, *out, "fasta", amino_acid_traits);
		n += seq.size();
	}
	out->close();
	delete out;
}