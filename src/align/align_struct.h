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

#ifndef ALIGN_STRUCT_H_
#define ALIGN_STRUCT_H_

#include "../basic/match.h"

struct local_match : public Hsp_data
{
	local_match() :
		query_anchor_(0),
		subject_(0)
	{ }
	local_match(int score) :
		Hsp_data(score)
	{ }
	local_match(int query_anchor, int subject_anchor, const Letter *subject, unsigned total_subject_len = 0) :
		total_subject_len_(total_subject_len),
		query_anchor_(query_anchor),
		subject_anchor(subject_anchor),
		subject_(subject)
	{ }
	local_match(unsigned len, unsigned query_begin, unsigned query_len, unsigned subject_len, unsigned gap_openings, unsigned identities, unsigned mismatches, signed subject_begin, signed score) :
		Hsp_data(score),
		query_anchor_(0),
		subject_(0)
	{ }
	unsigned total_subject_len_;
	signed query_anchor_, subject_anchor;
	const Letter *subject_;
};

#endif