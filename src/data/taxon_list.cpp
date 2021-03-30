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
#include "taxon_list.h"
#include "taxonomy.h"
#include "../util/log_stream.h"
#include "../util/io/text_input_file.h"
#include "../basic/config.h"
#include "../util/string/tokenizer.h"
#include "../util/algo/external_sort.h"
#include "../util/algo/sort_helper.h"

using std::set;
using std::endl;
using std::make_pair;

TaxonList::TaxonList(Deserializer &in, size_t size, size_t data_size):
	CompactArray<vector<uint32_t>>(in, size, data_size)
{}

static int mapping_file_format(const string& header) {
	string field1, field2;
	Util::String::Tokenizer tok(header, "\t");
	tok >> field1 >> field2;
	if (field1 == "accession" && field2 == "accession.version") {
		tok >> field1 >> field2;
		if (field1 == "taxid" && field2 == "gi" && !tok.good())
			return 0;
	}
	else if (field1 == "accession.version" && field2 == "taxid" && !tok.good())
		return 1;
	throw std::runtime_error("Accession mapping file header has to be in one of these formats:\naccession\taccession.version\ttaxid\tgi\naccession.version\ttaxid");
}

static void load_mapping_file(ExternalSorter<pair<string, uint32_t>>& sorter)
{
	unsigned taxid;
	TextInputFile f(config.prot_accession2taxid);
	f.getline();
	int format = mapping_file_format(f.line);
	string accession, last;

	while (!f.eof() && (f.getline(), !f.line.empty())) {
		if (format == 0)
			Util::String::Tokenizer(f.line, "\t") >> Util::String::Skip() >> accession >> taxid;
		else
			Util::String::Tokenizer(f.line, "\t") >> accession >> taxid;

		if (accession.empty())
			throw std::runtime_error("Empty accession field in line " + std::to_string(f.line_count));

		size_t i = accession.find(":PDB=");
		if (i != string::npos)
			accession.erase(i);

		i = accession.find_last_of('.');
		if (i != string::npos)
			accession.erase(i);

		if (accession != last)
			sorter.push(make_pair(accession, (uint32_t)taxid));

		last = accession;
	}
	f.close();
}

void TaxonList::build(OutputFile &db, ExternalSorter<pair<string, uint32_t>>& acc2oid, size_t seqs, Table& stats)
{
	typedef pair<string, uint32_t> T;

	task_timer timer("Loading taxonomy mapping file");
	ExternalSorter<T> acc2taxid;
	load_mapping_file(acc2taxid);

	timer.go("Joining accession mapping");
	acc2taxid.init_read();
	acc2oid.init_read();
	const auto cmp = [](const T& x, const T& y) { return x.first < y.first; };
	const auto value = [](const T& x, const T& y) { return make_pair(x.second, y.second); };
	auto it = join_sorted_lists(acc2oid, acc2taxid, cmp, value);

	ExternalSorter<pair<uint32_t, uint32_t>> oid2taxid;
	size_t acc_matched = 0;
	while (it.good()) {
		oid2taxid.push(*it);
		++it;
		++acc_matched;
	}

	timer.go("Writing taxon id list");
	db.set(Serializer::VARINT);
	oid2taxid.init_read();
	const uint32_t n = (uint32_t)seqs;
	auto taxid_it = merge_keys(oid2taxid, First<uint32_t, uint32_t>(), Second<uint32_t, uint32_t>(), 0);
	size_t mapped_seqs = 0;
	while (taxid_it.key() < n) {
		set<uint32_t> tax_ids = *taxid_it;
		tax_ids.erase(0);
		db << tax_ids;
		++taxid_it;
		if (!tax_ids.empty())
			++mapped_seqs;
	}
	timer.finish();

	stats("Accessions in database", acc2oid.count());
	stats("Entries in accession to taxid file", acc2taxid.count());
	stats("Database accessions mapped to taxid" , acc_matched);
	stats("Database sequences mapped to taxid", mapped_seqs);
}
