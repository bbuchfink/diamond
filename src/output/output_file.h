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

#ifndef OUTPUT_FILE_H_
#define OUTPUT_FILE_H_

#include <string>
#include "output.h"

using std::string;

struct Block_output : public Buffered_file
{

	struct Iterator
	{
		unsigned block_;
		bool same_subject_;
		Intermediate_record info_;
		bool operator<(const Iterator &rhs) const
		{ return info_.query_id > rhs.info_.query_id ||
				(info_.query_id == rhs.info_.query_id && (rhs.same_subject_ ||
						(!rhs.same_subject_ && info_.score < rhs.info_.score))); }
	};

	bool next(Iterator &it, unsigned subject, unsigned query)
	{
		if(this->eof())
			return false;
		it.info_.read(*this);
		it.block_ = block_;
		it.same_subject_ = it.info_.subject_id == subject && it.info_.query_id == query;
		return true;
	}

	Block_output(unsigned ref_block, const Temp_file &tmp_file):
		Buffered_file (tmp_file),
		block_ (ref_block)
	{ }

private:

	const unsigned block_;

};

#endif /* OUTPUT_FILE_H_ */
