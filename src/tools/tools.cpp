#include "tsv_record.h"
#include "../basic/config.h"
#include "../util/math/sparse_matrix.h"

using std::cout;
using std::endl;

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