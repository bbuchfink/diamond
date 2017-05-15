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

#ifndef MATCH_FUNC_H_
#define MATCH_FUNC_H_

#include "../basic/match.h"

inline int blast_frame(unsigned frame)
{ return frame <= 2 ? (int)frame+1 : 2-(int)frame; }

inline unsigned query_translated_begin(unsigned query_begin, unsigned frame, unsigned dna_len, bool query_translated)
{
	if(!query_translated)
		return query_begin;
	int f = frame <= 2 ? frame+1 : 2-frame;
	if (f > 0)
		return (query_begin - (f-1))/3;
	else
		return (dna_len + f - query_begin)/3;
}

#endif /* MATCH_FUNC_H_ */
