#include "reference.h"

String_set<char, '\0'>* ref_ids::data_ = nullptr;
Partitioned_histogram ref_hst;
unsigned current_ref_block;
SequenceSet* ref_seqs::data_ = nullptr;
SequenceSet* ref_seqs_unmasked::data_ = nullptr;
bool blocked_processing;
std::vector<uint32_t> block_to_database_id;
