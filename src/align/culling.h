#pragma once
#include "../basic/match.h"
#include "../output/output_format.h"

namespace Extension {

bool filter_hsp(Hsp& hsp,
	int source_query_len,
	const char *query_title,
	int subject_len,
	const char* subject_title,
	const Sequence& query_seq,
	const Sequence& subject_seq,
	const double query_self_aln_score,
	const double target_self_aln_score,
	const OutputFormat* output_format);

}