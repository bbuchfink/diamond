/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

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
