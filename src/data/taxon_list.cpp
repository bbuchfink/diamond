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

#include "taxon_list.h"
#include "taxonomy.h"

void TaxonList::build(OutputFile &db, FileBackedBuffer &accessions, size_t seqs)
{
	accessions.rewind();
	vector<string> a;
	Serializer out(db, Serializer::VARINT);
	vector<unsigned> t;
	for (size_t i = 0; i < seqs; ++i) {
		accessions >> a;
		t.clear();
		for (vector<string>::const_iterator j = a.begin(); j < a.end(); ++j)
			t.push_back(taxonomy.get(Taxonomy::Accession(j->c_str())));
		out << t;
	}
}