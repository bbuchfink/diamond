/****
DIAMOND protein aligner
Copyright (C) 2020 QIAGEN A/S (Aarhus, Denmark)
Code developed by Patrick Ettenhuber <patrick.ettenhuber@qiagen.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "output_format.h"
#include "recursive_parser.h"
#include "../data/queries.h"

void Clustering_format::print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const 
{
}

void Clustering_format::print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out) 
{
	out.write((uint32_t) query_block_to_database_id[r.query_id]);
	out.write((uint32_t) r.orig_subject_id);
	RecursiveParser rp(&r, format->c_str(), false);
	out.write(rp.evaluate());
}
