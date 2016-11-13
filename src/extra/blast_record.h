#ifndef BLAST_RECORD_H_
#define BLAST_RECORD_H_

#include <stdint.h>

using std::string;

struct blast_match {

	string query;
	string subject;
	double expect, bitscore, id;
	bool hit, stop;
	unsigned len, rid, ungapped_len, n, raw_score;

	bool empty() const
	{
		return query.length() == 0;
	}

	void set_empty()
	{
		query.resize(0);
	}

	bool operator<(const blast_match &rhs) const
	{
		return subject < rhs.subject;
	}

};

#endif /* BLAST_RECORD_H_ */
