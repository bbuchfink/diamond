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

#include <tuple>
#include <set>
#include "taxon_list.h"
#include "taxonomy.h"
#include "../util/log_stream.h"
#include "../util/io/text_input_file.h"
#include "../basic/config.h"
#include "../util/string/tokenizer.h"

using std::set;
using std::endl;
using std::tuple;
using std::get;
using std::make_tuple;

Serializer& operator<<(Serializer& s, const tuple<string, uint32_t>& x) {
	return s;
}

Serializer& operator<<(Serializer& s, const tuple<uint32_t, uint32_t>& x) {
	return s;
}

Deserializer& operator>>(Deserializer& s, tuple<string, uint32_t> x) {
	return s;
}

#include "../util/algo/external_sort.h"

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

static void load_mapping_file(ExternalSorter<tuple<string, uint32_t>, 0>& sorter)
{
	unsigned taxid;
	TextInputFile f(config.prot_accession2taxid);
	f.getline();
	int format = mapping_file_format(f.line);
	string accession;

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

		sorter.push(make_tuple(accession, (uint32_t)taxid));
	}
	f.close();
}

void TaxonList::build(OutputFile &db, ExternalSorter<tuple<string, uint32_t>, 0>& acc2oid, size_t seqs)
{
	ExternalSorter<tuple<string, uint32_t>, 0> acc2taxid;
	load_mapping_file(acc2taxid);

	task_timer timer("Joining accession mapping");
	TupleJoinIterator<tuple<string, uint32_t>, tuple<string, uint32_t>, 0, 0> it(acc2oid, acc2taxid);
	ExternalSorter<tuple<uint32_t, uint32_t>, 0> oid2taxid;
	while (it.good()) {
		auto p = *it;
		oid2taxid.push(make_tuple(get<1>(p.first), get<1>(p.second)));
		++it;
	}

	
	/*vector<string> a;
	db.set(Serializer::VARINT);
	set<unsigned> t;
	size_t mapped = 0, mappings = 0, len_errors = 0;
	for (size_t i = 0; i < seqs; ++i) {
		accessions >> a;
		for (vector<string>::const_iterator j = a.begin(); j < a.end(); ++j) {
			try {
				t.insert(taxonomy.get(Taxonomy::Accession(j->c_str())));
			}
			catch (AccessionLengthError &) {
				++len_errors;
			}
		}
		t.erase(0);
		db << t;
		mappings += t.size();
		if (!t.empty())
			++mapped;
		t.clear();
	}
	timer.finish();
	message_stream << mapped << " sequences mapped to taxonomy, " << mappings << " total mappings." << endl;
	if (len_errors)
		message_stream << "Warning: " << len_errors << " sequences ignored due to accession length overflow." << endl;*/
}
