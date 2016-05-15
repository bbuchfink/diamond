#include "queries.h"

unsigned current_query_chunk;
Sequence_set* query_source_seqs::data_ = 0;
Sequence_set* query_seqs::data_ = 0;
String_set<0>* query_ids::data_ = 0;
auto_ptr<seed_histogram> query_hst;