#include "daa_record.h"

const string & DAA_match_record::query_name() const
{
	return parent_.query_name;
}

const vector<Letter>& DAA_match_record::query() const
{
	return parent_.context[frame];
}

uint64_t DAA_match_record::db_letters() const
{
	return parent_.file_.db_letters();
}

unsigned DAA_match_record::query_end() const
{
	if (parent_.file_.mode() == mode_blastp) {
		return query_begin + translated_query_len - 1;
	}
	else if (parent_.file_.mode() == mode_blastx) {
		int len = (int)translated_query_len * 3 * (frame>2 ? -1 : 1);
		return (int)query_begin + (len > 0 ? -1 : 1) + len;
	}
	else if (parent_.file_.mode() == mode_blastn) {
		int len = (int)translated_query_len*(frame>0 ? -1 : 1);
		return (int)query_begin + (len > 0 ? -1 : 1) + len;
	}
	else
		return 0;
}
