/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include <set>
#include "output_format.h"
#include "../data/taxonomy.h"
#include "../util/sequence/sequence.h"

using std::set;
using std::vector;

void Taxon_format::print_match(const HspContext &r, Output::Info &info)
{
	const vector<TaxId> taxons(info.db->taxids(r.subject_oid));
	if (taxons.empty())
		return;
	evalue = std::min(evalue, r.evalue());
	try {
		for (vector<TaxId>::const_iterator i = taxons.begin(); i != taxons.end(); ++i)
			taxid = info.db->taxon_nodes().get_lca(taxid, *i);
	}
	catch (std::exception &) {
		fprintf(stderr, "Query=%s Subject=%s\n", r.query_title.c_str(), r.target_title.c_str());
		throw;
	}
}

void Taxon_format::print_query_epilog(Output::Info &info) const
{
	info.out.write_until(info.query.title, Util::Seq::id_delimiters);
	info.out << '\t' << taxid << '\t';
	if (taxid != 0)
		info.out.print_e(evalue);
	else
		info.out << '0';
	info.out << '\n';
}