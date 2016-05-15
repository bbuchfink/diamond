#include "reference.h"

String_set<0>* ref_ids::data_ = 0;
Ref_map ref_map;
seed_histogram ref_hst;
unsigned current_ref_block;
Reference_header ref_header;
Sequence_set* ref_seqs::data_ = 0;
