/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include <set>
#include "taxon_list.h"
#include "taxonomy.h"
#include "../util/log_stream.h"

using std::set;

TaxonList::TaxonList(Deserializer &in, size_t size, size_t data_size):
	CompactArray<vector<unsigned> >(in, size, data_size)
{}

void TaxonList::build(OutputFile &db, FileBackedBuffer &accessions, size_t seqs)
{
	task_timer timer("Writing taxon id lists");
	vector<string> a;
	db.set(Serializer::VARINT);
	set<unsigned> t;
	size_t mapped = 0, mappings = 0;
	for (size_t i = 0; i < seqs; ++i) {
		accessions >> a;
		for (vector<string>::const_iterator j = a.begin(); j < a.end(); ++j)
			t.insert(taxonomy.get(Taxonomy::Accession(j->c_str())));
		t.erase(0);
		db << t;
		mappings += t.size();
		if (!t.empty())
			++mapped;
		t.clear();
	}
	timer.finish();
	message_stream << mapped << " sequences mapped to taxonomy, " << mappings << " total mappings." << endl;
}