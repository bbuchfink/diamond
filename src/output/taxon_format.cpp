/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include <set>
#include "output_format.h"
#include "util/sequence/sequence.h"
#include "data/sequence_file.h"

using std::set;
using std::vector;
using std::string;

static string taxon_lineage(TaxId taxid, SequenceFile& db) {
	const vector<TaxId> lin = db.lineage(taxid);
	if (lin.empty())
		return "N/A";
	string out(db.taxon_scientific_name(lin.front()));
	for (vector<TaxId>::const_iterator i = lin.begin() + 1; i != lin.end(); ++i)
		out += "; " + db.taxon_scientific_name(*i);
	return out;
}

TaxonFormat::TaxonFormat() :
	OutputFormat(taxon, HspValues::NONE, Output::Flags::DEFAULT_REPORT_UNALIGNED),
	taxid(0),
	evalue(std::numeric_limits<double>::max())
{
	needs_taxon_id_lists = true;
	needs_taxon_nodes = true;
	needs_taxon_scientific_names = config.include_lineage;
}

void TaxonFormat::print_match(const HspContext &r, Output::Info &info)
{
	const vector<TaxId> taxons(info.db->taxids(r.subject_oid));
	if (taxons.empty())
		return;
	evalue = std::min(evalue, r.evalue());
	for (vector<TaxId>::const_iterator i = taxons.begin(); i != taxons.end(); ++i)
		taxid = info.db->get_lca(taxid, *i);
}

void TaxonFormat::print_query_epilog(Output::Info &info) const
{
	info.out.write_until(info.query.title, Util::Seq::id_delimiters);
	info.out << '\t' << taxid << '\t';
	if (taxid > 0)
		info.out.print_e(evalue);
	else
		info.out << '0';
	if (config.include_lineage)
		info.out << '\t' << (taxid > 0 ? taxon_lineage(taxid, *info.db) : "N/A");
	info.out << '\n';
}