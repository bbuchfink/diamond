#include <vector>
#include "tsv_record.h"
#include "../basic/config.h"
#include "../util/math/sparse_matrix.h"
#include "../basic/value.h"
#include "../util/util.h"
#include "../basic/sequence.h"

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