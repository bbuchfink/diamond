/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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
#include "util/sequence/sequence.h"

void PAF_format::print_query_intro(Output::Info& info) const
{
	if (info.unaligned) {
		info.out.write_until(info.query.title, Util::Seq::id_delimiters);
		info.out << "\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
	}
}

void PAF_format::print_match(const HspContext& r, Output::Info& info)
{
	info.out.write_until(r.query_title.c_str(), Util::Seq::id_delimiters);
	info.out << '\t' << r.query.source().length() << '\t'
		<< r.query_source_range().begin_ << '\t'
		<< r.query_source_range().end_ - 1 << '\t'
		<< (Frame(r.frame()).strand == FORWARD ? '+' : '-') << '\t';
	print_title(info.out, r.target_title.c_str(), false, false, "<>");

	info.out << '\t' << r.subject_len << '\t'
		<< r.subject_range().begin_ << '\t'
		<< r.subject_range().end_ - 1 << '\t'
		<< r.identities() << '\t'
		<< r.length() << '\t';
#ifdef WITH_DNA
    if(config.command == ::Config::blastn){
        info.out << r.mapping_quality() << '\t'
        << "cm:i:" << r.n_anchors() << '\t';
    } else
        info.out << "255" << '\t';
#else
    info.out << "255" << '\t';
#endif
    info.out	<< "AS:i:" << (uint32_t)score_matrix.bitscore(r.score()) << '\t'
		<< "ZR:i:" << r.score() << '\t'
		<< "ZE:f:";
	info.out.print_e(r.evalue());
	info.out << '\n';
}