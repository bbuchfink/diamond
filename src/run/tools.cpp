#include "tools.h"
#include "../basic/config.h"
#include "../data/sequence_set.h"
#include "../util/seq_file_format.h"
#include "../data/queries.h"
#include "../data/load_seqs.h"

void get_seq()
{
	const Sequence_file_format *format_n(guess_format(config.query_file));
	Compressed_istream query_file(config.query_file);
	load_seqs(query_file, *format_n, &query_seqs::data_, query_ids::data_, query_source_seqs::data_, std::numeric_limits<size_t>::max());
	cout << query_ids::get()[config.seq_no].c_str() << endl;
}